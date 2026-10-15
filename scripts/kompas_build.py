#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Построение венца нутационной передачи в КОМПАС-3D v21 через COM API.

ВНИМАНИЕ: скрипт написан без доступа к КОМПАС и НЕ запускался на реальном API.
Проверена только логика Python на заглушках. Вызовы API (особенно операции по
сечениям, плоскость под углом, круговой массив) могут потребовать правки.
Поэтому построение разбито на этапы, есть режим диагностики и ключи ориентации.

Нужно: Windows, 64-разрядный Python, `pip install pywin32`, запущенный или
устанавливаемый КОМПАС-3D v21. Данные берутся из JSON, который делает
scripts/nutation_gears.py (models/z1_z2/kompas_*.json). Единицы -- мм.

Построение (система координат детали: вершина конусов O в начале, ось детали Z):
  1. заготовка: вращение контура (две сферы R_in/R_out и два конуса) вокруг оси Z;
  2. впадина: вырез по сечениям между двумя плоскостями, перпендикулярными образующей
     делительного конуса (все рёбра проходят через O, поэтому боковые стороны зуба --
     точные конические поверхности);
  3. круговой массив выреза вокруг оси Z на Z зубьев.

Использование:
  python kompas_build.py --diag                       # что есть в вашем КОМПАС
  python kompas_build.py models/z1_z2/kompas_z1_wheel.json --stage blank
  python kompas_build.py models/z1_z2/kompas_z1_wheel.json --stage gap
  python kompas_build.py models/z1_z2/kompas_z1_wheel.json --stage full --save z1_wheel.m3d
"""
import argparse
import json
import math
import re
import sys

# GUID библиотек типов. API7 и константы 2D подтверждены по открытому коду на Python
# для КОМПАС; API5 и константы 3D взяты из памяти -- при ошибке см. --diag.
API5_GUIDS = ["{0422828C-F174-495E-AC5D-D31014DBBE87}"]
CONST3D_GUIDS = ["{2CAF168C-7961-4B90-9DA2-701419BEEFE3}"]
CONST2D_GUIDS = ["{75C9F5D0-B5B8-4526-8681-9903C567D2ED}"]


# ------------------------------------------------------------------ подключение
def _typelib(guids):
    from win32com.client import gencache
    last = None
    for g in guids:
        try:
            mod = gencache.EnsureModule(g, 0, 1, 0)
            if mod is not None:
                return mod
        except Exception as e:                       # noqa: BLE001
            last = e
    raise RuntimeError("не удалось загрузить библиотеку типов %s: %s" % (guids, last))


def connect():
    import pythoncom
    from win32com.client import Dispatch
    mod5 = _typelib(API5_GUIDS)
    api5 = mod5.KompasObject(Dispatch("Kompas.Application.5")._oleobj_.QueryInterface(
        mod5.KompasObject.CLSID, pythoncom.IID_IDispatch))
    api5.Visible = True
    return api5, _typelib(CONST3D_GUIDS).constants


def diag():
    """Печатает найденные библиотеки типов КОМПАС и подходящие имена констант."""
    try:
        import winreg
        print("Библиотеки типов с 'kompas' в описании:")
        with winreg.OpenKey(winreg.HKEY_CLASSES_ROOT, "TypeLib") as root:
            for i in range(winreg.QueryInfoKey(root)[0]):
                guid = winreg.EnumKey(root, i)
                try:
                    with winreg.OpenKey(root, guid) as k:
                        for j in range(winreg.QueryInfoKey(k)[0]):
                            ver = winreg.EnumKey(k, j)
                            desc = winreg.QueryValue(k, ver)
                            if re.search("kompas|компас", desc, re.I):
                                print("  %s  v%s  %s" % (guid, ver, desc))
                except OSError:
                    pass
    except Exception as e:                           # noqa: BLE001
        print("реестр прочитать не удалось:", e)
    try:
        api5, c = connect()
        print("API5 подключён.")
        names = [n for n in dir(c) if n.startswith("o3d_")]
        pat = re.compile("loft|offset|angle|rotat|circular|axis|plane(XO|YO)|sketch|section", re.I)
        print("Константы o3d_* (%d), подходящие по теме:" % len(names))
        for n in names:
            if pat.search(n):
                print("  %s = %s" % (n, getattr(c, n)))
    except Exception as e:                           # noqa: BLE001
        print("подключение не удалось:", type(e).__name__, e)


def const(c, *names):
    """Первая существующая константа из списка имён."""
    for n in names:
        if hasattr(c, n):
            return getattr(c, n)
    close = [x for x in dir(c) if any(n.lower()[:8] in x.lower() for n in names)][:12]
    raise RuntimeError("нет константы %s; похожие: %s (см. --diag)" % (names, close))


# ------------------------------------------------------------------- примитивы
class Builder:
    def __init__(self, api5, c, opts):
        self.api5, self.c, self.o = api5, c, opts
        self.doc = self.part = None

    def new_part(self):
        self.doc = self.api5.Document3D()
        self.doc.Create(False, True)                 # видимый документ, деталь
        self.part = self.doc.GetPart(-1)             # -1 -- верхний компонент
        return self.part

    def default(self, name):
        return self.part.GetDefaultEntity(const(self.c, name))

    def offset_plane(self, base, h, label):
        e = self.part.NewEntity(const(self.c, "o3d_planeOffset"))
        d = e.GetDefinition()
        d.direction = True
        d.offset = h
        d.SetPlane(base)
        e.name = label
        self._create(e, label)
        return e

    def angle_plane(self, base, axis, deg, label):
        e = self.part.NewEntity(const(self.c, "o3d_planeAngle"))
        d = e.GetDefinition()
        d.angle = deg
        d.SetPlane(base)
        d.SetAxis(axis)
        e.name = label
        self._create(e, label)
        return e

    def sketch(self, plane, label, draw):
        e = self.part.NewEntity(const(self.c, "o3d_sketch"))
        d = e.GetDefinition()
        d.SetPlane(plane)
        e.name = label
        self._create(e, label)
        doc2d = d.BeginEdit()
        try:
            draw(doc2d)
        finally:
            d.EndEdit()
        return e

    @staticmethod
    def _create(entity, label):
        if not entity.Create():
            raise RuntimeError("Create() вернул False: %s" % label)

    # --- 2D
    @staticmethod
    def polyline(doc2d, pts, closed=True, style=1):
        n = len(pts)
        for i in range(n if closed else n - 1):
            (x1, y1), (x2, y2) = pts[i], pts[(i + 1) % n]
            doc2d.ksLineSeg(x1, y1, x2, y2, style)

    @staticmethod
    def arc3(doc2d, p1, pm, p2, style=1):
        doc2d.ksArcBy3Points(p1[0], p1[1], pm[0], pm[1], p2[0], p2[1], style)


# ------------------------------------------------------------------ построение
def build_blank(b, prof):
    """Заготовка: вращение контура между сферами R_in/R_out и конусами c_a, c_b."""
    r_in, r_out = prof["r_in"], prof["r_out"]
    ca, cb = (math.radians(prof["blank"][k]) for k in ("c_a_deg", "c_b_deg"))
    fz = -1.0 if b.o.flip_z else 1.0

    def pt(r, c):                                     # (радиус от оси, высота вдоль оси)
        return (r * math.sin(c), fz * r * math.cos(c))

    cm = 0.5 * (ca + cb)
    plane = b.default("o3d_planeXOZ")

    def draw(doc2d):
        b.polyline(doc2d, [pt(r_in, ca), pt(r_out, ca)], closed=False)           # конус c_a
        b.arc3(doc2d, pt(r_out, ca), pt(r_out, cm), pt(r_out, cb))               # внешняя сфера
        b.polyline(doc2d, [pt(r_out, cb), pt(r_in, cb)], closed=False)           # конус c_b
        b.arc3(doc2d, pt(r_in, cb), pt(r_in, cm), pt(r_in, ca))                  # внутренняя сфера
        doc2d.ksLineSeg(0.0, -1.5 * r_out, 0.0, 1.5 * r_out, 3)                  # ось вращения

    sk = b.sketch(plane, "Эскиз заготовки", draw)
    e = b.part.NewEntity(const(b.c, "o3d_bossRotated"))
    d = e.GetDefinition()
    d.SetSketch(sk)
    d.SetSideParam(True, 360.0)
    b._create(e, "Вращение заготовки")
    return e


def gap_planes(b, prof):
    """Две плоскости, перпендикулярные образующей: XOY, повёрнутая вокруг оси Y на delta,
    затем смещённая на h_a и h_b."""
    delta = prof["delta_deg"] * (-1.0 if b.o.tilt_sign < 0 else 1.0)
    tilted = b.angle_plane(b.default("o3d_planeXOY"), b.default("o3d_axisOY"), delta, "Плоскость под образующую")
    pa = b.offset_plane(tilted, prof["gap"]["h_a"], "Сечение A")
    pb = b.offset_plane(tilted, prof["gap"]["h_b"], "Сечение B")
    return pa, pb


def build_gap(b, prof):
    xy = prof["gap"]["unit_xy"]
    sx = -1.0 if b.o.flip_x else 1.0
    pa, pb = gap_planes(b, prof)
    sketches = []
    for plane, h, label in ((pa, prof["gap"]["h_a"], "Контур впадины A"),
                            (pb, prof["gap"]["h_b"], "Контур впадины B")):
        pts = [(sx * x * h, y * h) for x, y in xy]
        sketches.append(b.sketch(plane, label, lambda d, pts=pts: b.polyline(d, pts)))
    e = b.part.NewEntity(const(b.c, "o3d_cutLoft"))
    d = e.GetDefinition()
    coll = d.Sketchs()
    coll.Clear()
    for s in sketches:
        coll.Add(s)
    b._create(e, "Вырез впадины по сечениям")
    return e


def build_pattern(b, prof, cut):
    e = b.part.NewEntity(const(b.c, "o3d_circularCopy"))
    d = e.GetDefinition()
    d.SetAxis(b.default("o3d_axisOZ"))
    d.SetCopyParamAlongDir(prof["z"], 360.0 / prof["z"], True, False)
    d.GetOperationArray().Add(cut)
    b._create(e, "Массив впадин")
    return e


STAGES = ("blank", "gap", "full")


def build_part(api5, c, prof, opts):
    b = Builder(api5, c, opts)
    b.new_part()
    steps = [("заготовка", lambda: build_blank(b, prof))]
    if opts.stage in ("gap", "full"):
        steps.append(("вырез впадины", lambda: build_gap(b, prof)))
    result = {}
    for label, fn in steps:
        print("  - %s ..." % label, end=" ")
        try:
            result[label] = fn()
        except Exception as ex:                       # noqa: BLE001
            print("ОШИБКА")
            raise RuntimeError("этап '%s': %s: %s" % (label, type(ex).__name__, ex)) from ex
        print("ок")
    if opts.stage == "full":
        print("  - круговой массив ...", end=" ")
        try:
            build_pattern(b, prof, result["вырез впадины"])
        except Exception as ex:                       # noqa: BLE001
            print("ОШИБКА")
            raise RuntimeError("этап 'массив': %s: %s" % (type(ex).__name__, ex)) from ex
        print("ок")
    if opts.save:
        b.doc.SaveAs(opts.save)
        print("сохранено:", opts.save)
    return b


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("profile", nargs="?", help="kompas_*.json из nutation_gears.py")
    ap.add_argument("--stage", choices=STAGES, default="blank")
    ap.add_argument("--diag", action="store_true")
    ap.add_argument("--save", help="путь для сохранения .m3d")
    ap.add_argument("--tilt-sign", type=int, choices=(-1, 1), default=1,
                    help="-1, если плоскости сечений наклонились не в ту сторону")
    ap.add_argument("--flip-x", action="store_true", help="зеркально по локальной оси x эскиза сечений")
    ap.add_argument("--flip-z", action="store_true", help="поменять знак высоты в эскизе заготовки")
    opts = ap.parse_args(argv)

    if opts.diag:
        diag()
        return 0
    if not opts.profile:
        ap.error("нужен файл профиля или --diag")
    with open(opts.profile, encoding="utf-8") as f:
        prof = json.load(f)
    print("Деталь %s: z=%d, делительный конус %.4f°, этап %s" %
          (prof["name"], prof["z"], prof["delta_deg"], opts.stage))
    api5, c = connect()
    build_part(api5, c, prof, opts)
    print("Готово. Проверьте дерево построения и модель в КОМПАС.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

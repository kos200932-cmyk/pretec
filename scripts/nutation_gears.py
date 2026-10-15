#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Нутационная (прецессирующая) передача с двухвенцовым сателлитом.

Схема: колесо Z1 -- венец Z2 сателлита, венец Z3 сателлита -- колесо Z4.
Передаточное число i = 1 / (1 - Z1*Z3 / (Z2*Z4)).

Модель зацепления
  * Оси колеса и сателлита пересекаются в точке O (начало координат) под углом нутации theta.
  * Делительные конусы пары имеют общую вершину O и касаются по образующей g.
  * Колёса (Z1, Z4) -- наружные зубья, венцы сателлита (Z2, Z3) -- внутренние зубья.
  * Зубья прямые, профиль -- сферическая эвольвента (точное сопряжение для пары
    с общей вершиной). Боковая поверхность зуба -- конус с вершиной O.
  * Торцы венцов -- сферы радиусов R_in и R_out (R_out -- внешний делительный
    конусный размер). Модуль m задан на внешнем торце.

Все длины в мм, углы во внутренних расчётах в радианах. Нужен только numpy.
"""
import argparse
import json
import os
import struct

import numpy as np


# ---------------------------------------------------------------- кинематика
def gear_ratio(z):
    z1, z2, z3, z4 = z
    return 1.0 / (1.0 - (z1 * z3) / (z2 * z4))


def pair_cones(z_w, z_s, theta):
    """Полуугол делительного конуса колеса (наружные зубья) и венца сателлита
    (внутренние зубья, z_s > z_w) при угле между осями theta."""
    u = z_s / z_w
    if u <= np.cos(theta):
        raise ValueError("Z_сат/Z_кол слишком мало для такого угла нутации")
    d_w = np.arctan2(np.sin(theta), u - np.cos(theta))
    return d_w, d_w + theta


# ------------------------------------------------------------ сферическая геометрия
def sph(az, c):
    """Единичные векторы по азимуту az и полярному углу (колатитюде) c."""
    return np.stack([np.sin(c) * np.cos(az), np.sin(c) * np.sin(az), np.cos(c)], axis=-1)


def to_sph(p):
    return np.arctan2(p[..., 1], p[..., 0]), np.arccos(np.clip(p[..., 2], -1.0, 1.0))


def rot_z(a):
    ca, sa = np.cos(a), np.sin(a)
    return np.array([[ca, -sa, 0.0], [sa, ca, 0.0], [0.0, 0.0, 1.0]])


def sat_frame(theta):
    """Столбцы -- оси x', y', z' системы сателлита в мировой системе.
    z' -- ось сателлита (наклонена на theta от оси Z), x' -- в сторону образующей g."""
    return np.array([[np.cos(theta), 0.0, -np.sin(theta)],
                     [0.0, 1.0, 0.0],
                     [np.sin(theta), 0.0, np.cos(theta)]])


def inv_az(gb, c):
    """Азимутальное смещение сферической эвольвенты с основным конусом gb
    на колатитюде c (c >= gb); в точке c = gb смещение равно 0."""
    cos_e = np.clip(np.cos(c) / np.cos(gb), -1.0, 1.0)
    e = np.arccos(cos_e)
    return e / np.sin(gb) - np.arctan(np.tan(e) / np.sin(gb))


# ------------------------------------------------------------------ профиль венца
def make_outline(z, delta, alpha, ha_ang, hf_ang, internal,
                 n_flank=32, n_tip=5, n_root=8, relief=0.04):
    """Контур венца в координатах (азимут, колатитюд) на сфере.
    Зуб с номером 0 расположен симметрично относительно азимута 0.
    Возвращает азимуты строго по возрастанию на полный оборот."""
    gb = np.arcsin(np.sin(delta) * np.cos(alpha))     # основной конус
    a_p = inv_az(gb, delta)
    s_half = np.pi / (2.0 * z)                         # половина толщины зуба на делительном конусе
    pitch = 2.0 * np.pi / z

    if internal:
        c_tip, c_root = delta - ha_ang, delta + hf_ang
        if c_tip < gb:
            raise ValueError("вершина зуба лежит внутри основного конуса")
        c_seq = np.linspace(c_root, c_tip, n_flank)    # от впадины к вершине
        sign = -1.0
    else:
        c_tip, c_root = delta + ha_ang, delta - hf_ang
        c_seq = np.linspace(max(c_root, gb), c_tip, n_flank)
        sign = 1.0

    az_l = -s_half + sign * (inv_az(gb, c_seq) - a_p)  # левая сторона зуба
    c_l = c_seq.copy()
    if (not internal) and c_root < gb:                 # участок ниже основного конуса
        az_l = np.insert(az_l, 0, az_l[0] - relief * pitch)
        c_l = np.insert(c_l, 0, c_root)

    tip = np.linspace(az_l[-1], -az_l[-1], n_tip + 2)[1:-1]
    az_r, c_r = -az_l[::-1], c_l[::-1]
    root = np.linspace(-az_l[0], az_l[0] + pitch, n_root + 2)[1:-1]

    az_t = np.concatenate([az_l, tip, az_r, root])
    c_t = np.concatenate([c_l, np.full(n_tip, c_tip), c_r, np.full(n_root, c_root)])
    az = np.concatenate([az_t + k * pitch for k in range(z)])
    c = np.tile(c_t, z)
    if not np.all(np.diff(az) > 0):
        raise ValueError("контур зуба неоднозначен по азимуту (подрез?)")
    # впадина между зубом 0 и зубом 1, центрированная на азимуте 0 (от вершины до вершины)
    gap_az = np.concatenate([az_r, root, az_l + pitch]) - pitch / 2.0
    gap_c = np.concatenate([c_r, np.full(n_root, c_root), c_l])
    return dict(az=az, c=c, gb=gb, a_p=a_p, c_tip=c_tip, c_root=c_root, internal=internal,
                gap_az=gap_az, gap_c=gap_c, pitch=pitch)


# ----------------------------------------------------------------- тело венца
def ring_solid(az, c_a, c_b, r_in, r_out, a_is_low):
    """Замкнутое тело между сферами r_in, r_out, ограниченное по колатитюду
    кривыми A (контур зубьев) и B (окружность ступицы/обода)."""
    n = len(az)
    i0 = np.arange(n)
    i1 = (i0 + 1) % n
    V = np.vstack([r_in * sph(az, c_a), r_in * sph(az, c_b),
                   r_out * sph(az, c_a), r_out * sph(az, c_b)])

    def ix(k, curve, i):
        return (2 * k + curve) * n + i

    low, high = (0, 1) if a_is_low else (1, 0)
    F = []
    a, b, c, d = ix(1, low, i0), ix(1, high, i0), ix(1, high, i1), ix(1, low, i1)   # внешний торец
    F += [np.stack([a, b, c], 1), np.stack([a, c, d], 1)]
    a, b, c, d = ix(0, low, i0), ix(0, high, i0), ix(0, high, i1), ix(0, low, i1)   # внутренний торец
    F += [np.stack([a, c, b], 1), np.stack([a, d, c], 1)]
    q00, q10, q11, q01 = ix(0, low, i0), ix(1, low, i0), ix(1, low, i1), ix(0, low, i1)
    F += [np.stack([q00, q10, q11], 1), np.stack([q00, q11, q01], 1)]               # нижняя стенка
    q00, q10, q11, q01 = ix(0, high, i0), ix(1, high, i0), ix(1, high, i1), ix(0, high, i1)
    F += [np.stack([q00, q11, q10], 1), np.stack([q00, q01, q11], 1)]               # верхняя стенка
    return V, np.vstack(F)


def mesh_stats(V, F):
    tri = V[F]
    vol = np.einsum('ij,ij->i', tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0
    e = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]).astype(np.int64)
    key = e[:, 0] * len(V) + e[:, 1]
    rkey = e[:, 1] * len(V) + e[:, 0]
    uniq, counts = np.unique(key, return_counts=True)
    closed = bool(counts.max() == 1 and np.isin(rkey, uniq).all())
    area = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    return dict(triangles=int(len(F)), vertices=int(len(V)), closed_manifold=closed,
                volume_mm3=float(vol), degenerate=int((area < 1e-12).sum()))


def write_stl(path, V, F, name="pretec"):
    tri = V[F]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = n / np.where(ln == 0, 1.0, ln)
    rec = np.zeros(len(F), dtype=np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')]))
    rec['n'], rec['v'] = n, tri
    with open(path, 'wb') as f:
        f.write(name.encode('ascii', 'ignore')[:80].ljust(80, b' '))
        f.write(struct.pack('<I', len(F)))
        f.write(rec.tobytes())


# ------------------------------------------------------- проверка зацепления
def _bound(az_q, az, c):
    a0 = az[0]
    return np.interp((az_q - a0) % (2 * np.pi) + a0, az, c)


def check_meshing(z_w, z_s, theta, W, S, c_rim, c_bore, r_out, n_steps=90):
    """Вращаем колесо на phi и сателлит на psi = phi*z_w/z_s вокруг их осей и ищем
    точки контура одного венца, попавшие в материал другого (в микрометрах)."""
    M = sat_frame(theta)
    psi0 = np.pi / z_s
    pw = sph(W['az'], W['c'])
    ps = sph(S['az'], S['c'])
    pen_w = pen_s = 0.0
    for phi in np.linspace(0.0, 2 * np.pi / z_w, n_steps, endpoint=False):
        a_s = psi0 + phi * z_w / z_s
        # точки колеса -> система сателлита: материал сателлита c >= граница
        loc = (pw @ rot_z(phi).T) @ M @ rot_z(a_s)
        az_l, c_l = to_sph(loc)
        d = c_l - _bound(az_l, S['az'], S['c'])
        d = np.where(c_l <= c_rim, d, -1.0)
        pen_w = max(pen_w, float(d.max()))
        # точки сателлита -> система колеса: материал колеса c <= граница
        wld = (ps @ rot_z(a_s).T @ M.T) @ rot_z(phi)
        az_l, c_l = to_sph(wld)
        d = _bound(az_l, W['az'], W['c']) - c_l
        d = np.where(c_l >= c_bore, d, -1.0)
        pen_s = max(pen_s, float(d.max()))
    return max(pen_w, 0.0) * r_out * 1e3, max(pen_s, 0.0) * r_out * 1e3


# ------------------------------------------------- данные для скрипта КОМПАС-3D
def gap_section(P, delta, c_open_shift):
    """Контур впадины (замкнутый) в плоскости, перпендикулярной образующей
    делительного конуса, на расстоянии 1 от вершины O вдоль образующей.
    Образующая лежит в плоскости XZ: g = (sin d, 0, cos d). Локальные оси плоскости:
    e1 = (cos d, 0, -sin d), e2 = (0, 1, 0). Сечение на расстоянии h -- это то же
    самое, умноженное на h (все ребра проходят через O)."""
    az, c = P['gap_az'], P['gap_c']
    c_open = P['c_tip'] + c_open_shift
    az = np.concatenate([az, [az[-1], az[0]]])             # замыкание со стороны открытого края
    c = np.concatenate([c, [c_open, c_open]])
    u = sph(az, c)
    g = np.array([np.sin(delta), 0.0, np.cos(delta)])
    e1 = np.array([np.cos(delta), 0.0, -np.sin(delta)])
    e2 = np.array([0.0, 1.0, 0.0])
    k = u @ g
    if k.min() < 0.3:
        raise ValueError("контур впадины слишком далёк от образующей для плоского сечения")
    p = u / k[:, None]
    return np.stack([p @ e1, p @ e2], axis=1)


def kompas_part(name, kind, z, delta, P, r_in, r_out, c_a, c_b):
    shift = 0.08 if kind == 'external' else -0.08
    xy = gap_section(P, delta, shift)
    return dict(
        name=name, kind=kind, z=int(z), delta_deg=float(np.degrees(delta)),
        pitch_deg=float(np.degrees(P['pitch'])),
        gap_center_az_deg=float(np.degrees(P['pitch'] / 2.0)),   # азимут центра впадины в STL
        r_in=float(r_in), r_out=float(r_out),
        blank=dict(c_a_deg=float(np.degrees(c_a)), c_b_deg=float(np.degrees(c_b))),
        gap=dict(h_a=float(0.8 * r_in), h_b=float(1.15 * r_out), unit_xy=xy.round(9).tolist()),
    )


def check_gap_section(P, delta, part):
    """Обратное преобразование точек контура в (азимут, колатитюд): ошибка относительно
    исходного контура и проверка, что ребра между сечениями проходят через O."""
    xy = np.array(part['gap']['unit_xy'])[:len(P['gap_az'])]
    g = np.array([np.sin(delta), 0.0, np.cos(delta)])
    e1 = np.array([np.cos(delta), 0.0, -np.sin(delta)])
    e2 = np.array([0.0, 1.0, 0.0])
    p = g + np.outer(xy[:, 0], e1) + np.outer(xy[:, 1], e2)
    u = p / np.linalg.norm(p, axis=1, keepdims=True)
    az, c = to_sph(u)
    err = max(np.abs(az - P['gap_az']).max(), np.abs(c - P['gap_c']).max())
    pa, pb = part['gap']['h_a'] * p, part['gap']['h_b'] * p
    cross = np.linalg.norm(np.cross(pa, pb), axis=1).max()   # 0, если точки коллинеарны с O
    return float(err), float(cross)


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--z', type=int, nargs=4, default=[15, 16, 18, 17], metavar=('Z1', 'Z2', 'Z3', 'Z4'))
    ap.add_argument('--theta', type=float, default=10.511009, help='угол нутации, градусы')
    ap.add_argument('--module', type=float, default=2.0, help='модуль на внешнем торце, мм')
    ap.add_argument('--alpha', type=float, default=20.0, help='угол профиля, градусы')
    ap.add_argument('--ha', type=float, default=1.0, help='коэффициент высоты головки')
    ap.add_argument('--hf', type=float, default=1.25, help='коэффициент высоты ножки')
    ap.add_argument('--width', type=float, default=1 / 3, help='ширина венца как доля внешнего конусного размера')
    ap.add_argument('--pair', choices=['12', '34'], default='12', help='какая пара строится')
    ap.add_argument('--out', default='models/z1_z2')
    ap.add_argument('--kompas-pts', type=int, default=10, help='точек на эвольвенте для скрипта КОМПАС')
    args = ap.parse_args()

    z = args.z
    theta = np.radians(args.theta)
    alpha = np.radians(args.alpha)
    print("Передаточное число i = %.6f" % gear_ratio(z))
    for nm, (zw, zs) in {'1-2': (z[0], z[1]), '3-4': (z[3], z[2])}.items():
        dw, ds = pair_cones(zw, zs, theta)
        print("  пара %s: δ колеса %.4f°, δ венца сателлита %.4f°" % (nm, np.degrees(dw), np.degrees(ds)))

    z_w, z_s = (z[0], z[1]) if args.pair == '12' else (z[3], z[2])
    names = ('z%d_wheel' % (1 if args.pair == '12' else 4), 'z%d_crown' % (2 if args.pair == '12' else 3))
    d_w, d_s = pair_cones(z_w, z_s, theta)
    R = args.module * z_w / (2 * np.sin(d_w))
    assert abs(R - args.module * z_s / (2 * np.sin(d_s))) < 1e-9, "конусные размеры пары не совпали"
    b = args.width * R
    r_in, r_out = R - b, R
    ha_ang, hf_ang = args.ha * args.module / R, args.hf * args.module / R

    W = make_outline(z_w, d_w, alpha, ha_ang, hf_ang, internal=False)
    S = make_outline(z_s, d_s, alpha, ha_ang, hf_ang, internal=True)

    c_bore = float(np.arcsin(0.4 * np.sin(d_w)))
    c_rim = S['c_root'] + 0.10
    Vw, Fw = ring_solid(W['az'], W['c'], np.full_like(W['az'], c_bore), r_in, r_out, a_is_low=False)
    Vs, Fs = ring_solid(S['az'], S['c'], np.full_like(S['az'], c_rim), r_in, r_out, a_is_low=True)

    stats = {names[0]: mesh_stats(Vw, Fw), names[1]: mesh_stats(Vs, Fs)}
    pen_w, pen_s = check_meshing(z_w, z_s, theta, W, S, c_rim, c_bore, r_out)

    # сборка в мировой системе: ось колеса -- Z, ось сателлита наклонена на theta
    M = sat_frame(theta)
    Vs_world = (Vs @ rot_z(np.pi / z_s).T) @ M.T
    Va = np.vstack([Vw, Vs_world])
    Fa = np.vstack([Fw, Fs + len(Vw)])

    os.makedirs(args.out, exist_ok=True)
    write_stl(os.path.join(args.out, names[0] + '.stl'), Vw, Fw, names[0])
    write_stl(os.path.join(args.out, names[1] + '.stl'), Vs, Fs, names[1])
    write_stl(os.path.join(args.out, 'pair_assembly.stl'), Va, Fa, 'pair_assembly')

    summary = dict(
        z=z, pair=args.pair, i=gear_ratio(z), theta_deg=args.theta, module_mm=args.module,
        alpha_deg=args.alpha, ha_coef=args.ha, hf_coef=args.hf,
        pitch_cone_wheel_deg=float(np.degrees(d_w)), pitch_cone_crown_deg=float(np.degrees(d_s)),
        base_cone_wheel_deg=float(np.degrees(W['gb'])), base_cone_crown_deg=float(np.degrees(S['gb'])),
        cone_distance_mm=float(R), face_width_mm=float(b),
        pitch_diameter_wheel_mm=args.module * z_w, pitch_diameter_crown_mm=args.module * z_s,
        mesh=stats, penetration_um=dict(wheel_into_crown=pen_w, crown_into_wheel=pen_s),
        files=[n + '.stl' for n in names] + ['pair_assembly.stl'],
    )
    with open(os.path.join(args.out, 'summary.json'), 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    # профили для скрипта КОМПАС-3D (облегчённый контур, kompas-pts точек на боковую сторону)
    Wk = make_outline(z_w, d_w, alpha, ha_ang, hf_ang, internal=False, n_flank=args.kompas_pts, n_tip=1, n_root=3)
    Sk = make_outline(z_s, d_s, alpha, ha_ang, hf_ang, internal=True, n_flank=args.kompas_pts, n_tip=1, n_root=3)
    parts = [kompas_part(names[0], 'external', z_w, d_w, Wk, r_in, r_out, c_bore, W['c_tip']),
             kompas_part(names[1], 'internal', z_s, d_s, Sk, r_in, r_out, S['c_tip'], c_rim)]
    for part, P, d in zip(parts, (Wk, Sk), (d_w, d_s)):
        part['gap_check'] = dict(zip(('max_angle_error_rad', 'apex_collinearity_mm2'), check_gap_section(P, d, part)))
        with open(os.path.join(args.out, 'kompas_%s.json' % part['name']), 'w', encoding='utf-8') as f:
            json.dump(part, f, ensure_ascii=False)
        xy = np.array(part['gap']['unit_xy'])
        print("  КОМПАС %s: %d точек контура впадины, сечение h=R: x %.2f..%.2f, y %.2f..%.2f мм, ошибка %.1e рад, коллинеарность %.1e" %
              (part['name'], len(xy), xy[:, 0].min() * r_out, xy[:, 0].max() * r_out,
               xy[:, 1].min() * r_out, xy[:, 1].max() * r_out,
               part['gap_check']['max_angle_error_rad'], part['gap_check']['apex_collinearity_mm2']))

    print("Пара %s: R = %.3f мм, ширина %.3f мм" % (args.pair, R, b))
    print("  основные конусы: колесо %.3f°, венец %.3f°" % (np.degrees(W['gb']), np.degrees(S['gb'])))
    for k, v in stats.items():
        print("  %s: %d треуг., замкнута=%s, объём %.1f мм³, вырожденных %d" %
              (k, v['triangles'], v['closed_manifold'], v['volume_mm3'], v['degenerate']))
    print("  проникновение профилей за оборот пары: колесо→венец %.2f мкм, венец→колесо %.2f мкм" % (pen_w, pen_s))


if __name__ == '__main__':
    main()

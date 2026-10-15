#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Геометрические данные для построения пары в КОМПАС-3D (без обращения к КОМПАС).

Как строится деталь (колесо или венец сателлита) в КОМПАС:
  1. Заготовка -- вращение замкнутого контура на плоскости XOZ вокруг оси Z.
     Контур ограничен двумя сферами (R_in, R_out) с центром в вершине O конусов
     и двумя коническими поверхностями.
  2. Впадина -- операция «вырезать по сечениям» между двумя эскизами на плоскостях
     z = h1 и z = h2. Боковая поверхность зуба -- конус с вершиной O, поэтому
     сечение на z = h2 это сечение на z = h1, умноженное на h2/h1 по x и y;
     линейчатая поверхность между ними точно совпадает с боковой поверхностью зуба.
  3. Остальные впадины -- то же сечение, повёрнутое вокруг оси Z на k * шаг.

Скрипт считает эти данные, проверяет их независимо (объём детали из данных сечений
сравнивается с объёмом сетки из nutation_gears.py) и пишет JSON.
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nutation_gears as ng  # noqa: E402

N_TIP, N_ROOT = 5, 8
_trapz = getattr(np, 'trapezoid', None) or np.trapz


def plane_pts(az, c, h):
    """Точки сферы (азимут, колатитюд) -> сечение плоскостью z = h вдоль лучей из O."""
    t = np.tan(c) * h
    return np.stack([t * np.cos(az), t * np.sin(az)], axis=-1)


def gap_indices(n_total, z, n_flank_total=None):
    period = n_total // z
    n_l = (period - N_TIP - N_ROOT) // 2
    return np.arange(n_l + N_TIP, period + n_l), n_l, period


def rotate(poly, a):
    ca, sa = math.cos(a), math.sin(a)
    p = np.asarray(poly)
    return np.stack([p[:, 0] * ca - p[:, 1] * sa, p[:, 0] * sa + p[:, 1] * ca], axis=1)


# ----------------------------------------------------------------- детали
def _part(name, z, outline, internal, r_in, r_out, c_a, c_b, ext=0.06, n_ext=8):
    az, c = outline['az'], outline['c']
    idx, n_l, period = gap_indices(len(az), z)
    az_g, c_g = az[idx], c[idx]
    c_tip = outline['c_tip']
    c_ext = c_tip - ext if internal else c_tip + ext
    az_a, az_b = az_g[0], az_g[-1]
    path_n = len(idx)
    clos_az = np.concatenate([[az_b], np.linspace(az_b, az_a, n_ext + 2)[1:-1], [az_a]])
    clos_c = np.full(len(clos_az), c_ext)
    paz = np.concatenate([az_g, clos_az])
    pc = np.concatenate([c_g, clos_c])
    z_lo = r_in * math.cos(pc.max())
    z_hi = r_out * math.cos(pc.min())
    h1, h2 = 0.85 * z_lo, 1.15 * z_hi
    poly1, poly2 = plane_pts(paz, pc, h1), plane_pts(paz, pc, h2)

    P = lambda r, cc: (r * math.sin(cc), r * math.cos(cc))
    d = math.degrees
    profile = [
        ["arc", 0.0, 0.0, r_in, d(math.pi / 2 - c_b), d(math.pi / 2 - c_a)],
        ["line", *P(r_in, c_a), *P(r_out, c_a)],
        ["arc", 0.0, 0.0, r_out, d(math.pi / 2 - c_b), d(math.pi / 2 - c_a)],
        ["line", *P(r_in, c_b), *P(r_out, c_b)],
    ]
    return dict(name=name, z=z, internal=internal, pitch_deg=360.0 / z,
                r_in=r_in, r_out=r_out, c_a=c_a, c_b=c_b, c_tip=c_tip,
                profile=profile, axis=[0.0, 0.0, 0.0, 1.2 * r_out],
                h1=h1, h2=h2, path_n=path_n,
                poly1=poly1.tolist(), poly2=poly2.tolist())


def pair_data(z=(15, 16, 18, 17), theta_deg=10.511009, module=2.0, alpha_deg=20.0,
              ha=1.0, hf=1.25, width=1 / 3, pair='12', n_flank=14):
    theta, alpha = math.radians(theta_deg), math.radians(alpha_deg)
    z_w, z_s = (z[0], z[1]) if pair == '12' else (z[3], z[2])
    nm_w, nm_s = ('z1_wheel', 'z2_crown') if pair == '12' else ('z4_wheel', 'z3_crown')
    d_w, d_s = ng.pair_cones(z_w, z_s, theta)
    R = module * z_w / (2 * math.sin(d_w))
    b = width * R
    r_in, r_out = R - b, R
    ha_ang, hf_ang = ha * module / R, hf * module / R

    def outlines(nf):
        return (ng.make_outline(z_w, d_w, alpha, ha_ang, hf_ang, False, n_flank=nf, n_tip=N_TIP, n_root=N_ROOT),
                ng.make_outline(z_s, d_s, alpha, ha_ang, hf_ang, True, n_flank=nf, n_tip=N_TIP, n_root=N_ROOT))

    W, S = outlines(n_flank)
    c_bore = float(np.arcsin(0.4 * np.sin(d_w)))
    c_rim = S['c_root'] + 0.10
    wheel = _part(nm_w, z_w, W, False, r_in, r_out, c_bore, W['c_tip'])
    crown = _part(nm_s, z_s, S, True, r_in, r_out, S['c_tip'], c_rim)

    # опорные объёмы из сеток того же генератора (n_flank = 32), независимый путь расчёта
    W32, S32 = outlines(32)
    Vw, Fw = ng.ring_solid(W32['az'], W32['c'], np.full_like(W32['az'], c_bore), r_in, r_out, a_is_low=False)
    Vs, Fs = ng.ring_solid(S32['az'], S32['c'], np.full_like(S32['az'], c_rim), r_in, r_out, a_is_low=True)
    wheel['ref_volume_mm3'] = ng.mesh_stats(Vw, Fw)['volume_mm3']
    crown['ref_volume_mm3'] = ng.mesh_stats(Vs, Fs)['volume_mm3']
    wheel['_dense'], crown['_dense'] = W32, S32
    meta = dict(z=list(z), pair=pair, theta_deg=theta_deg, module=module, alpha_deg=alpha_deg,
                ha=ha, hf=hf, width=width, n_flank=n_flank, cone_distance_mm=R, face_width_mm=b)
    return meta, wheel, crown


# --------------------------------------------------------------- проверки
def _seg_cross(p, q, r, s):
    o = lambda a, b, c: (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return o(p, q, r) * o(p, q, s) < 0 and o(r, s, p) * o(r, s, q) < 0


def is_simple(poly):
    n = len(poly)
    for i in range(n):
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            if _seg_cross(poly[i], poly[(i + 1) % n], poly[j], poly[(j + 1) % n]):
                return False
    return True


def _dist_to_polyline(pts, poly):
    a, b = poly[:-1], poly[1:]
    ab = b - a
    t = np.clip(((pts[:, None, :] - a) * ab).sum(-1) / (ab * ab).sum(-1), 0.0, 1.0)
    proj = a + t[..., None] * ab
    return np.linalg.norm(pts[:, None, :] - proj, axis=-1).min(axis=1)


def verify(pd):
    """Независимые проверки данных одной детали. Возвращает словарь результатов."""
    res = {}
    p1, p2 = np.array(pd['poly1']), np.array(pd['poly2'])
    res['sections_similar'] = bool(np.allclose(p2, p1 * (pd['h2'] / pd['h1']), atol=1e-9))
    res['polygon_simple'] = is_simple(p1.tolist()) and is_simple(p2.tolist())
    seg = np.linalg.norm(np.roll(p2, -1, axis=0) - p2, axis=1)
    res['min_segment_mm'] = float(seg.min())
    res['n_segments'] = int(len(p1))

    # объём детали из данных сечений: заготовка минус z впадин
    path = p2[:pd['path_n']]
    dense = np.concatenate([np.linspace(path[i], path[i + 1], 60, endpoint=False) for i in range(len(path) - 1)]
                           + [path[-1:]])
    p3 = np.column_stack([dense, np.full(len(dense), pd['h2'])])
    u = p3 / np.linalg.norm(p3, axis=1, keepdims=True)
    az, c = np.arctan2(u[:, 1], u[:, 0]), np.arccos(u[:, 2])
    omega_gap = float(_trapz(np.abs(np.cos(c) - math.cos(pd['c_tip'])), az))
    omega_blank = 2 * math.pi * abs(math.cos(pd['c_a']) - math.cos(pd['c_b']))
    shell = (pd['r_out'] ** 3 - pd['r_in'] ** 3) / 3.0
    vol = shell * (omega_blank - pd['z'] * omega_gap)
    res['volume_from_sections_mm3'] = vol
    res['volume_mesh_mm3'] = pd['ref_volume_mm3']
    res['volume_rel_diff'] = abs(vol - pd['ref_volume_mm3']) / pd['ref_volume_mm3']

    # отклонение хорд сечения от «истинной» кривой (плотный контур, 32 точки на сторону)
    dn = pd['_dense']
    idx, _, _ = gap_indices(len(dn['az']), pd['z'])
    true_pts = plane_pts(dn['az'][idx], dn['c'][idx], pd['h2'])
    res['chord_dev_um'] = float(_dist_to_polyline(true_pts, path).max() * 1e3)
    return res


def clean(pd):
    return {k: v for k, v in pd.items() if not k.startswith('_')}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--z', type=int, nargs=4, default=[15, 16, 18, 17])
    ap.add_argument('--theta', type=float, default=10.511009)
    ap.add_argument('--module', type=float, default=2.0)
    ap.add_argument('--alpha', type=float, default=20.0)
    ap.add_argument('--ha', type=float, default=1.0)
    ap.add_argument('--hf', type=float, default=1.25)
    ap.add_argument('--width', type=float, default=1 / 3)
    ap.add_argument('--pair', choices=['12', '34'], default='12')
    ap.add_argument('--n-flank', type=int, default=14, help='точек на боковую сторону зуба в эскизе')
    ap.add_argument('--out', default='models/z1_z2/kompas_data.json')
    a = ap.parse_args()

    meta, wheel, crown = pair_data(tuple(a.z), a.theta, a.module, a.alpha, a.ha, a.hf, a.width, a.pair, a.n_flank)
    for pd in (wheel, crown):
        r = verify(pd)
        pd['checks'] = r
        print("%s: z=%d, сечения подобны=%s, контур без самопересечений=%s, сегментов %d, мин. сегмент %.4f мм"
              % (pd['name'], pd['z'], r['sections_similar'], r['polygon_simple'], r['n_segments'], r['min_segment_mm']))
        print("   плоскости h1=%.3f, h2=%.3f мм; отклонение хорд от кривой %.2f мкм"
              % (pd['h1'], pd['h2'], r['chord_dev_um']))
        print("   объём: из сечений %.2f мм³, сетка %.2f мм³, расхождение %.4f%%"
              % (r['volume_from_sections_mm3'], r['volume_mesh_mm3'], 100 * r['volume_rel_diff']))
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w', encoding='utf-8') as f:
        json.dump(dict(meta=meta, parts=[clean(wheel), clean(crown)]), f, ensure_ascii=False, indent=1)
    print("записано:", a.out)


if __name__ == '__main__':
    main()

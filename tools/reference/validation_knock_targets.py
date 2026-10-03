"""Validation round 1 (knock & octane rating): experimental targets for the knock-limited CR test.

Reads only committed CFR fixtures (tools/reference/cfr_data_*.py outputs) and writes
test/fixtures/validation_knock_targets.json:

  * ASTM guide-table compression ratios (rigid-raise counter -> CR) for PRF 40-100, RON and MON (9/16 in
    venturi), plus the least-squares slope dCR/dON over 60-100 (D2699-15a / D2700-14 Tables A4.1).
  * The measured STANDARD-KNOCK STATE at RON (what a physically faithful model must reproduce at the
    guide-table CR, independent of any KI proxy): combustion phasing and knock point
    (Kolodziej & Wallner 2017 Fig. 9), peak pressure (Choi 2018 Fig. 9 PRF100 -13 deg flush-mount),
    knock-point pressure, post-KP pressure-rise rate, MAPO and KP position in the heat release
    (Rockstroh 2018 standard point), MAPO (Hoth 2021 Tables 6-7), knock-point scatter (Pal 2018 Fig. 9).
  * Measured trapped fresh charge from fuel flow and lambda (KW17 Fig. 2, Choi 2018 Table 6) -> the
    volumetric efficiency the model's end-gas density/temperature must be consistent with.
  * Spark-advance equivalence of one octane number near standard knock (Pal 2018 Fig. 1 KU slope,
    ASTM D2699 A2.5.2 spread 12-15 KU/ON).
  * Intake-temperature insensitivity of the knock-limited CR of PRF90 at 1.0 bar (Rockstroh 2018 Fig. 5).
  * The 501-C detonation-meter INPUT transfer function: least-squares fit of V(t) = K*HP_fc[LP_6.5kHz[p]]
    to Hoth 2021a Fig. 2 (pressure + meter input of the same cycle): the meter input is the pressure
    signal high-passed at fc ~ 50 Hz (i.e. the D-1 dp/dt signal integrated by a ~50 Hz RC stage), so
    the meter responds to the band-limited pressure development after the knock point, not to the
    6/10/14 kHz ringing (consistent with Swarts via Rockstroh 2018 p. 2 and p. 8).

Run: .venv/bin/python tools/reference/validation_knock_targets.py
"""
import json
import os

import numpy as np
from scipy import optimize, signal

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
FIX = os.path.join(ROOT, 'test', 'fixtures')


def load(name):
    with open(os.path.join(FIX, name)) as f:
        return json.load(f)


def ds(d, key):
    return next(x for x in d['datasets'] if x['key'] == key)


astm = load('cfr_astm_guide_tables.json')
ONS = [40, 50, 60, 70, 75, 80, 85, 90, 93, 95, 97, 98, 100]


def guide(method, on):
    t = astm['RON_9_16']['onToCounter'] if method == 'RON' else astm['MON']['9/16']['onToCounter']
    i = t['octaneNumber'].index(float(on))
    return {'counter': t['counter'][i], 'compressionRatio': t['compressionRatioRigidRaise'][i]}


guideTable = {m: {str(on): guide(m, on) for on in ONS} for m in ('RON', 'MON')}
slopes = {}
for m in ('RON', 'MON'):
    x = np.array([60, 70, 80, 85, 90, 95, 100], float)
    y = np.array([guideTable[m][str(int(o))]['compressionRatio'] for o in x])
    slopes[m] = float(np.polyfit(x, y, 1)[0])

# ---- standard-knock state (RON) ----
kw = load('cfr_kw2017_ron98.json')
kwm = ds(kw, 'kw2017_ron98_operating_metrics')
fig = {f['figure']: f['values'] for f in kwm['figures']}
choi = load('cfr_choi2018_traces.json')
c13 = ds(choi, 'choi2018_fig9_prf100_st_m13_flushMount')['derivedFromTrace']
rock = load('cfr_rockstroh2018_knock_vs_cr.json')['datasets'][0]
rsp = rock['standardKnockPoint']
hoth = load('cfr_hoth_knock_metrics.json')
h6 = ds(hoth, 'hoth2021a_fig6_ku_mapo_vs_lambda')['reported']
pal = load('cfr_pal2018_knock.json')
p9 = ds(pal, 'pal2018_fig9_knock_points')
p1 = ds(pal, 'pal2018_fig1_ku_mapo_vs_spark')['data']

standardKnockState = {
    'KW2017_PRF98': {
        'ca10Deg': fig['Fig. 9']['CA10']['PRF98'],
        'ca50Deg': fig['Fig. 9']['CA50']['PRF98'],
        'knockPointDeg': fig['Fig. 9']['knockPoint']['PRF98'],
        'peakPressureRiseRatePaPerDeg': fig['Fig. 11']['PPRR_PaPerDeg']['PRF98'],
        'knockPressurePeakPa': fig['Fig. 11']['KP_PK_Pa']['PRF98'],
        'gIMEPPa': fig['Fig. 4']['gIMEP']['PRF98'],
        'lambda': fig['Fig. 2']['lambda']['PRF98'],
        'source': 'Kolodziej & Wallner 2017 Figs. 2, 4, 9, 11 (OSTI 1394801); spark-plug transducer',
    },
    'Choi2018_PRF100_m13_flush': {
        'peakPressurePa': c13['peakPressurePa'],
        'peakPressureDeg': c13['peakPressureThetaDeg'],
        'compressionRatio': 7.55,
        'source': 'Choi et al. 2018 Fig. 9 (OSTI 1501884), 300-cycle average',
    },
    'Rockstroh2018_PRF90_standard': {
        'compressionRatio': rsp['compressionRatio'],
        'pressureAtKnockPointPa': rsp['digitisedFromFigures']['pressureAtKnockPointPa']['value'],
        'pressureRiseRateAfterKPPaPerDeg': rsp['digitisedFromFigures']['pressureRiseRateAfterKPPaPerDeg']['value'],
        'peakFilteredPressurePa': rsp['digitisedFromFigures']['peakFilteredPressurePa']['value'],
        'mapoPa': rsp['digitisedFromFigures']['mapoPa']['value'],
        'knockPointFractionOfHeatRelease': 0.6,
        'source': 'Rockstroh et al. 2018 (SAE 2018-01-0210) standard point, Figs. 2, 5, 6; KP at ~60 % of aHR',
    },
    'Hoth2021_MAPO_4to18kHz': {
        'PRF98_standardRON_sparkPlug': h6['table6_PRF98_standardRON']['sparkPlugTransducer_18kHz']['mapoBar'] * 1e5,
        'PRF98_standardRON_knockmeterPort': h6['table6_PRF98_standardRON']['kistler6044A_in_knockmeter_port_18kHz']['mapoBar'] * 1e5,
        'PRF95_at_RON95_CR_52.8KU': h6['table7']['PRF95']['peakMapoBar'] * 1e5,
        'unit': 'Pa',
        'source': 'Hoth & Kolodziej 2021 Part 1 Tables 6-7 (OSTI 1880351)',
    },
    'Pal2018_knockPoint': {'meanDeg': 11.03, 'sdDeg': 0.97, 'source': 'Pal et al. 2018 Fig. 9 (219 cycles)'},
    'endGasMassFractionAtKnockPoint': {
        'value': [0.3, 0.4],
        'note': 'KP at ~60-70 % MFB at standard RON conditions (Rockstroh 2018 Fig. 2 and conclusions: '
                'intake heating moves KP from ~70 to 40 % MFB at 50 KU); KW17: KP 3.2 deg after CA50.',
    },
}

# ---- measured trapped fresh charge ----
AFST = 15.1  # KW17 Table 1 (PRF98)
fr = fig['Fig. 2']['fuelRate']['PRF98']
lam = fig['Fig. 2']['lambda']['PRF98']
cycles_per_s = 600 / 60 / 2
fuel_cycle = fr / cycles_per_s
air_cycle = fuel_cycle * AFST * lam
# ambient dry air density at 52 C (RON IAT), 1 atm
rho_dry = 101325 / (287.05 * 325.15)
Vd = np.pi / 4 * 0.08255 ** 2 * 0.1143
tab6 = ds(choi, 'choi2018_tables')['table6_PRF98_standardRON_default']
measuredCharge = {
    'KW2017_PRF98': {'fuelPerCycleKg': fuel_cycle, 'airPerCycleKg': air_cycle,
                     'volumetricEfficiencyDryAirAmbient52C': air_cycle / (rho_dry * Vd)},
    'Choi2018_PRF98_Table6': {'fuelPerCycleKg': tab6['fuelRateKgPerH'] / 3600 / cycles_per_s,
                              'airPerCycleKg': tab6['fuelRateKgPerH'] / 3600 / cycles_per_s * AFST * tab6['lambda']},
    'Pal2018_TPA_trappedMassKg': 0.000628,
    'Choi2018_TPA_TIVC_K_regression': '0.5385*T_port[C] + 394.34 (~410 K at T_port 28 C)',
    'Choi2018_TPA_Tu_minus20_K_regression': '0.6567*T_port[C] + 678.42 (~697 K at T_port 28 C)',
    'note': 'AFst 15.1 (KW17 Table 1); fresh charge from measured fuel flow and lambda; the TPA values are model outputs.',
}

# ---- spark-advance equivalence of one ON ----
sa = np.array(p1['sparkAdvanceDegBTDC'])
ku = np.array(p1['knockUnits'])
sel = (sa > 10.5) & (sa < 14.5)
kuPerDeg = float(np.polyfit(sa[sel], ku[sel], 1)[0])
sparkEquivalence = {
    'kuPerDegAdvance_Pal2018_Fig1_10.8to14deg': kuPerDeg,
    'kuPerOn_ASTM_spread': [12, 15],
    'octaneNumbersPerDegAdvance': [kuPerDeg / 15, kuPerDeg / 12],
}

# ---- Rockstroh: knock-limited CR vs intake temperature (PRF90, lambda 1, 1.0 bar) ----
def cr_at_50(series):
    k = series['knockUnits']
    cr = np.array(k['compressionRatio'])
    v = np.array(k['value'])
    o = np.argsort(cr)
    cr, v = cr[o], v[o]
    for i in range(1, len(cr)):
        if v[i - 1] < 50 <= v[i]:
            return float(cr[i - 1] + (50 - v[i - 1]) / (v[i] - v[i - 1]) * (cr[i] - cr[i - 1]))
    return None


rdat = rock['data']
rockstrohIntakeT = {}
for name in ('Tin33C_Pin1.00bar', 'Tin90C_Pin1.00bar', 'Tin150C_Pin1.00bar', 'Tin33C_Pin1.28bar', 'Tin150C_Pin1.28bar'):
    s = rdat[name]
    # Tin90C/1.0 bar: the lone 49.7 KU marker at CR 5.32 is a digitisation mis-assignment (KU 1.5-1.8
    # at CR 5.69-6.52); use the monotone part above CR 5.6.
    if name == 'Tin90C_Pin1.00bar':
        k = s['knockUnits']
        keep = [i for i, c in enumerate(k['compressionRatio']) if c > 5.6]
        s = {'knockUnits': {'compressionRatio': [k['compressionRatio'][i] for i in keep], 'value': [k['value'][i] for i in keep]}}
    c = cr_at_50(s)
    if c is None and name == 'Tin90C_Pin1.00bar':
        c = 6.93  # 44.4 KU at CR 6.90; extrapolated with the 33 C / 150 C slopes (~90 KU per CR)
    rockstrohIntakeT[name] = {'intakeMixtureTemperatureK': rdat[name]['intakeMixtureTemperatureK'],
                              'intakePressurePa': rdat[name]['intakePressurePa'], 'compressionRatioAt50KU': c}

# ---- 501-C input transfer function (Hoth 2021a Fig. 2) ----
h2 = ds(hoth, 'hoth2021a_fig2_prf100_detonation_meter_input')['data']
th = np.array(h2['thetaDeg'])
p = np.array(h2['pPa']) / 1e5
g = np.arange(-359.9, 359.9, 0.1)
P = np.interp(g, th, p)
V = np.interp(g, np.array(h2['meterInputThetaDeg']), np.array(h2['meterInputV']))
dt = 0.1 / 3600.0  # 600 rpm
fs = 1 / dt
P3 = np.tile(P, 3)
bl, al = signal.butter(1, 6500 / (fs / 2), 'low')


def model(par):
    K, fc, off = par
    b, a = signal.butter(1, fc / (fs / 2), 'high')
    y = signal.lfilter(b, a, P3)[len(P):2 * len(P)]
    return K * signal.lfilter(bl, al, y) + off


fit = optimize.least_squares(lambda q: model(q) - V, [0.02, 50, 0], bounds=([0, 0.1, -1], [1, 5000, 1]))
knockmeterInput = {
    'model': 'V = K * HP1_fc[LP1_6500Hz[p]] + offset (first-order Butterworth, periodic cycle)',
    'K_V_per_bar': float(fit.x[0]), 'highPassCornerHz': float(fit.x[1]), 'offsetV': float(fit.x[2]),
    'rmsResidualV': float(np.sqrt(np.mean(fit.fun ** 2))), 'peakV': float(V.max()),
    'source': 'Hoth & Kolodziej 2021 Part 1 Fig. 2 (OSTI 1880351), PRF100 knocking cycle',
}

out = {
    'description': 'Knock / octane-rating validation targets (validation round 1). Generated by '
                   'tools/reference/validation_knock_targets.py from committed cfr_* fixtures.',
    'generator': 'tools/reference/validation_knock_targets.py',
    'guideTable': guideTable,
    'guideTableSlopeCRperON_60to100': slopes,
    'standardKnockState': standardKnockState,
    'measuredCharge': measuredCharge,
    'sparkEquivalence': sparkEquivalence,
    'rockstrohIntakeTemperature': rockstrohIntakeT,
    'knockmeterInput': knockmeterInput,
}
with open(os.path.join(FIX, 'validation_knock_targets.json'), 'w') as f:
    json.dump(out, f, indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ('guideTable',)}, indent=1))

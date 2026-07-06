import numpy as np
d = 7.48e-3          # espesor de la placa [m]
lambda_ = 532e-9     # longitud de onda [m]

# Lista de datos del experimento: [m, phi_i en grados]
datos = [
    [121, 9.0],
    [142, 10.0],
    [160, 11.0],
    [181, 12.0],
    [202, 12.5],
    [222, 13.0],
    [243, 13.5],
    [266, 14.0],
    [286, 14.5],
    [307, 15.0],
    [330, 15.5],
]
def calcular_m(n, phi_i_deg, d, lambda_):
    phi_i = np.radians(phi_i_deg)           # ángulo de incidencia en radianes
    
    # Ley de Snell: sin(phi_i) = n * sin(phi_r)
    sin_phi_r = np.sin(phi_i) / n
    if abs(sin_phi_r) > 1:
        return np.nan                       
    phi_r = np.arcsin(sin_phi_r)            # ángulo de refracción en radianes
    
    # Ecuación (1) de Fendley
    term1 = n * (1/np.cos(phi_r) - 1)
    term2 = 1 - np.cos(phi_i - phi_r)/np.cos(phi_r)
    m_calc = (2*d / lambda_) * (term1 + term2)
    
    return m_calc, phi_r

def encontrar_n(m_exp, phi_i_deg, d, lambda_, n_min=1.3, n_max=1.6, paso=1e-5):
    n_opt = None
    m_opt = None
    phi_r_opt = None
    mejor_diferencia = np.inf
    
    n_prueba = n_min
    while n_prueba <= n_max:
        m_calc, phi_r = calcular_m(n_prueba, phi_i_deg, d, lambda_)
        if not np.isnan(m_calc):
            diff = abs(m_calc - m_exp)
            if diff < mejor_diferencia:
                mejor_diferencia = diff
                n_opt = n_prueba
                m_opt = m_calc
                phi_r_opt = phi_r
        n_prueba += paso
    
    return n_opt, phi_r_opt, m_opt

print("="*70)
print(f"d = {d*1000:.2f} mm,  lambda = {lambda_*1e9:.1f} nm")
print("="*70)
print(f"{'phi_i (°)':>10} | {'m (exp)':>8} | {'phi_r (°)':>10} | {'n':>10}")
print("-"*70)

resultados = []
for m_exp, phi_i_deg in datos:
    n, phi_r, m_calc = encontrar_n(m_exp, phi_i_deg, d, lambda_)
    phi_r_deg = np.degrees(phi_r) if phi_r is not None else np.nan
    resultados.append([phi_i_deg, m_exp, phi_r_deg, n])
    print(f"{phi_i_deg:10.1f} | {m_exp:8.0f} | {phi_r_deg:10.2f} | {n:10.6f}")

print("="*70)

n_values = [r[3] for r in resultados if r[3] is not None]
n_prom = np.mean(n_values)
n_std = np.std(n_values, ddof=1)
n_error = n_std / np.sqrt(len(n_values))

print(f"\nPromedio n = {n_prom:.6f} ± {n_error:.6f}")
print(f"Desviación estándar = {n_std:.6f}")
print(f"Número de datos = {len(n_values)}")
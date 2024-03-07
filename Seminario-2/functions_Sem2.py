from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
from IPython.display import SVG
import py3Dmol
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from scipy.stats import linregress

def MolTo3DView(smiles, size=(300, 300), style="stick", surface=False, opacity=0.5):
    """Draw molecule in 3D

    Args:
    ----
        mol: rdMol, molecule to show
        size: tuple(int, int), canvas size
        style: str, type of drawing molecule
               style can be 'line', 'stick', 'sphere', 'carton'
        surface, bool, display SAS
        opacity, float, opacity of surface, range 0.0-1.0
    Return:
    ----
        viewer: py3Dmol.view, a class for constructing embedded 3Dmol.js views in ipython notebooks.
    """
    assert style in ('line', 'stick', 'sphere', 'carton')
    mol = Chem.MolFromSmiles(smiles)
    mol_opt = optimize_conf(mol)
    mblock = Chem.MolToMolBlock(mol_opt)
    viewer = py3Dmol.view(width=size[0], height=size[1])
    viewer.addModel(mblock, 'mol')
    viewer.setStyle({style:{}})
    if surface:
        viewer.addSurface(py3Dmol.SAS, {'opacity': opacity})
    viewer.zoomTo()
    return viewer

def optimize_conf(mol):
    if mol is not None:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        return mol
    else:
        return None

def construct_molecule_2d(smiles):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)

    return mol

def construct_molecule_2d(smiles):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    return mol

def hydrolysis_products(smiles):
    # Convertir SMILES a objeto molécula
    mol = Chem.MolFromSmiles(smiles)

    # Inicializar lista para almacenar productos de hidrólisis
    hydrolysis_products_list = []

    # Verificar si la molécula tiene un grupo éster
    ester_pattern = Chem.MolFromSmarts('[CX3:1](=[O:2])[OX2:3][C:4]')
    amide_pattern = Chem.MolFromSmarts('[CX3:1](=[O:2])[NX3,NX2,#7]')#[C:4]')
    #amide_pattern = Chem.MolFromSmarts('[CX3:1](=[O:2])[NX3;H1:3]')#[C:4]')

    if mol.HasSubstructMatch(ester_pattern):
        if not mol.HasSubstructMatch(amide_pattern):
            print('Posee un éster!')
            # Definir reacción de hidrólisis para éster
            rxn_smarts = '[CX3:1](=[O:2])[OX2:3][C:4].[OX2;H2:5]>>[CX3:1](=[O:2])[OX2;H1:3].[C:4][OX2;H1:5]'
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            products = rxn.RunReactants((mol, Chem.MolFromSmiles('O')))
            for product in products:
                for p in product:
                    product_smiles = Chem.MolToSmiles(p)
                    hydrolysis_products_list.append(product_smiles)

    # Verificar si la molécula tiene un grupo amida
    if mol.HasSubstructMatch(amide_pattern) and not mol.HasSubstructMatch(ester_pattern):
      print('Posee una amida!')

      # Primera expresión de rxn_smarts
      rxn_smarts_1 = '[CX3:1](=[O:2])[NX3,NX2,#7;H1:3].[OX2;H2:5]>>[CX3:1](=[O:2])[OX2;H1:5].[NX3,NX2,#7;H2:3]'

      # Segunda expresión de rxn_smarts
      rxn_smarts_2 = '[CX3:1](=[O:2])[NX3,NX2,#7:3].[OX2;H2:5]>>[CX3:1](=[O:2])[OX2;H1:5].[NX3,NX2,#7;H1:3]'

      # Intentar con la primera expresión
      try:
          rxn = AllChem.ReactionFromSmarts(rxn_smarts_1)
          products = rxn.RunReactants((mol, Chem.MolFromSmiles('O')))

          if not products:
              # Si products está vacío, intentar con la segunda expresión
              rxn = AllChem.ReactionFromSmarts(rxn_smarts_2)
              products = rxn.RunReactants((mol, Chem.MolFromSmiles('O')))

              if not products:
                  print("Ambas expresiones de rxn_smarts no produjeron productos.")
              else:
                  #print("Se utilizó la segunda expresión exitosamente.")
                  # Resto del código utilizando rxn_smarts_2
                  for product in products:
                      for p in product:
                          product_smiles = Chem.MolToSmiles(p)
                          hydrolysis_products_list.append(product_smiles)
          else:
              #print("Se utilizó la primera expresión exitosamente.")
              # Resto del código utilizando rxn_smarts_1
              for product in products:
                  for p in product:
                      product_smiles = Chem.MolToSmiles(p)
                      hydrolysis_products_list.append(product_smiles)

      except Exception as e:
          print(f"Error al intentar la expresión: {e}")
          print("No se pudo construir rxn_smarts")



    if mol.HasSubstructMatch(ester_pattern):
        if mol.HasSubstructMatch(amide_pattern):
            print('Posee tanto éster como amida!')
            # Definir reacción de hidrólisis para éster
            rxn_smarts = '[CX3:1](=[O:2])[OX2:3][C:4].[OX2;H2:5]>>[CX3:1](=[O:2])[OX2;H1:3].[C:4][OX2;H1:5]'
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            products = rxn.RunReactants((mol, Chem.MolFromSmiles('O')))
            for product in products:
                for p in product:
                    product_smiles = Chem.MolToSmiles(p)
                    hydrolysis_products_list.append(product_smiles)
            for product_smiles in hydrolysis_products_list:
                product_mol = Chem.MolFromSmiles(product_smiles)
                rxn_smarts_amide = '[CX3:1](=[O:2])[NX3;H1:3].[OX2;H2:5]>>[CX3:1](=[O:2])[OX2;H1:5].[NX3;H2:3]'
                rxn_amide = AllChem.ReactionFromSmarts(rxn_smarts_amide)
                products_amide = rxn_amide.RunReactants((product_mol, Chem.MolFromSmiles('O')))
                for product_amide in products_amide:
                    for p_amide in product_amide:
                        product_amide_smiles = Chem.MolToSmiles(p_amide)
                        hydrolysis_products_list.append(product_amide_smiles)

    print("Reactivo")
    mol2d = construct_molecule_2d(smiles)
    display(mol2d)
    viewer = MolTo3DView(smiles, surface=False)
    viewer.show()
    print("Producto(s)")
    for unique_product in hydrolysis_products_list:
        mol2d = construct_molecule_2d(unique_product)
        display(mol2d)
        viewer = MolTo3DView(unique_product, surface=False)
        viewer.show()

    #return hydrolysis_products_list

def oxidize_molecule(input_smiles):
    mol = Chem.MolFromSmiles(input_smiles)

    # Verificar si la molécula tiene un grupo éster
    catecol_pattern1 = Chem.MolFromSmarts('[CX3,#6:1]([OH:2])=[#6:3]([OH:4])')
    catecol_pattern2 = Chem.MolFromSmarts('[CX3,#6:1]([OH:2])[#6:3]([OH:4])')

    # Inicializar conjuntos de productos
    unique_products1 = set()
    unique_products2 = set()

    if mol.HasSubstructMatch(catecol_pattern1) or mol.HasSubstructMatch(catecol_pattern2):
        #print('es catecol')

        # Intentar con la primera reacción SMARTS
        rxn_smarts1 = '[CX3,#6:1]([OH:2])=[CX3,#6:3]([OH:4])>>[CX3,#6:1](=[O:2])[CX3,#6:3](=[O:4])'
        rxn1 = AllChem.ReactionFromSmarts(rxn_smarts1)
        products1 = rxn1.RunReactants((mol,))
        for product1 in products1:
            for p1 in product1:
                product_smiles1 = Chem.MolToSmiles(p1)
                #print(product_smiles1)
                unique_products1.add(product_smiles1)

        # Verificar si la primera reacción no generó productos y, en ese caso, intentar con la segunda reacción
        if not unique_products1:
            rxn_smarts2 = '[CX3,#6:1]([OH:2])[CX3,#6:3]([OH:4])>>[CX3,#6:1](=[O:2])[CX3,#6:3](=[O:4])'
            rxn2 = AllChem.ReactionFromSmarts(rxn_smarts2)
            products2 = rxn2.RunReactants((mol,))
            for product2 in products2:
                for p2 in product2:
                    product_smiles2 = Chem.MolToSmiles(p2)
                    #print(product_smiles2)
                    unique_products2.add(product_smiles2)
    # Mostrar las imágenes 2D y 3D de los productos únicos
    for unique_product in unique_products1.union(unique_products2):
        print("Reactivo")
        mol2d = construct_molecule_2d(input_smiles)
        display(mol2d)
        viewer = MolTo3DView(input_smiles, surface=False)
        viewer.show()
        print("Producto(s)")
        mol2d = construct_molecule_2d(unique_product)
        display(mol2d)
        viewer = MolTo3DView(unique_product, surface=False)
        viewer.show()

def preparar_dataframe(arrhenius_data, temp_unit='C'):
    # Crear DataFrame original
    arrhenius_data_df = pd.DataFrame(arrhenius_data)

    # Convertir temperatura a Kelvin si es necesario
    if temp_unit == 'C':
        arrhenius_data_df['Temperature_K'] = arrhenius_data_df['Temperature'] + 273.15
    elif temp_unit == 'K':
        arrhenius_data_df['Temperature_K'] = arrhenius_data_df['Temperature']

    # Agregar columnas al nuevo DataFrame
    arrhenius_data_df['ln_k'] = np.log(arrhenius_data_df['k'])
    arrhenius_data_df['1_over_T'] = 1 / arrhenius_data_df['Temperature_K']

    return arrhenius_data_df

def arrhenius_equation(T, A, Ea_over_R):
    return A * np.exp(-Ea_over_R / T)

def linealizar_arrhenius(T, k):
    # Linealizar ln(k) = -Ea/R * 1/T + ln(A)
    inv_T = 1 / T
    ln_k = np.log(k)

    # Ajuste lineal
    slope, intercept, r_value, p_value, std_err = linregress(inv_T, ln_k)

    # Pendiente representa -Ea/R
    Ea_over_R = -slope
    Ea = Ea_over_R*1.987

    return slope, intercept, Ea_over_R, r_value, Ea

def graficar_arrhenius(data_df):
    sns.set(style="white")  # Configuración del estilo de Seaborn

    # Linealizar los datos
    slope, intercept, Ea_over_R, r_value, Ea = linealizar_arrhenius(data_df['Temperature_K'].values, data_df['k'].values)

    # Graficar los datos y ajuste lineal
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=data_df['1_over_T'], y=np.log(data_df['k']), label='Datos experimentales', color='blue')

    # Corregir la generación de la línea de ajuste lineal
    arrhenius_line = slope * data_df['1_over_T'].values + intercept
    plt.plot(data_df['1_over_T'].values, arrhenius_line, color='red', label='Ajuste lineal')

    plt.xlabel('1 / Temperatura (1/K)')
    plt.ylabel('ln(k)')
    plt.legend()

    plt.title(f'Ajuste Lineal de la Ecuación de Arrhenius\n$E_a/R$: {Ea_over_R:.2f}, $E_a$: {Ea:.2f}, $R^2$: {r_value**2:.2f}')

    # Imprimir valores de Ea y ordenada al origen
    print(f'Ea/R: {Ea_over_R:.4f}')
    print(f'Ea: {Ea:.4f}')
    print(f'Ordenada al origen (ln(A)): {intercept:.4f}')
    # Calcular e^x
    resultado = np.exp(intercept)

    # Imprimir el resultado
    print(f'A = {resultado}')

    plt.show()

def calcular_Ea(T1, T2, k1, k2, R=1.987):
    return R * np.log(k2 / k1) / (1 / T1 - 1 / T2)

def calcular_k(Ea, T1, T2, k1, R=1.987):
    return k1 * np.exp((Ea/R) * ((1 / T1) - (1 / T2)))

def calcular_concentracion(absorbancia, coef_extincion=1, longitud_camino=1):
    concentracion = absorbancia / (longitud_camino * coef_extincion)
    return concentracion

def generar_dataframe_desde_Abs(df, coef_extincion=1, longitud_camino=1):
    if 'concentracion' not in df.columns:
        df['concentracion'] = calcular_concentracion(df['absorbancia'], coef_extincion, longitud_camino)
        df['ln_concentracion'] = np.log(df['concentracion'])
        df['inversa_concentracion'] = 1 / df['concentracion']
    print(df)
    return df

def cinetica_orden_cero(t, c0, k):
    return c0 - k * t

def ajustar_cinetica_orden_cero(tiempo, concentracion):
    parametros_iniciales = [concentracion.iloc[0], 0.001]
    parametros_optimizados, matriz_covarianza = curve_fit(cinetica_orden_cero, tiempo, concentracion, p0=parametros_iniciales)

    concentracion_ajustada = cinetica_orden_cero(tiempo, *parametros_optimizados)
    r2 = r2_score(concentracion, concentracion_ajustada)
    sse = np.sum(np.square(concentracion - concentracion_ajustada))

    return parametros_optimizados[0], parametros_optimizados[1], r2, sse

def cinetica_orden_uno(t, c0, k):
    return c0 - k * t

def ajustar_cinetica_orden_uno(tiempo, concentracion_ln):
    slope, intercept, r_value, p_value, std_err = linregress(tiempo, concentracion_ln)
    r2 = r_value**2
    sse = np.sum(np.square(concentracion_ln - (slope * tiempo + intercept)))

    return intercept, -slope, r2, sse

def cinetica_orden_dos(t, c0, k):
    return (1 / c0) + k * t

def ajustar_cinetica_orden_dos(tiempo, inversa_concentracion):
    parametros_iniciales = [1 / inversa_concentracion.iloc[0], 0.001]
    parametros_optimizados, matriz_covarianza = curve_fit(cinetica_orden_dos, tiempo, inversa_concentracion, p0=parametros_iniciales)

    concentracion_ajustada = cinetica_orden_dos(tiempo, *parametros_optimizados)
    r2 = r2_score(inversa_concentracion, concentracion_ajustada)
    sse = np.sum(np.square(inversa_concentracion - concentracion_ajustada))

    return parametros_optimizados[0], parametros_optimizados[1], r2, sse

def graficar_cineticas(tiempo, datos, parametros_optimizados, orden, color_puntos='blue', color_linea='red'):
    tiempo = np.array(tiempo)
    datos = np.array(datos)

    plt.subplot(1, 3, orden)
    sns.scatterplot(x=tiempo, y=datos, label='Datos experimentales', color=color_puntos)
    plt.xlabel('Tiempo')

    if orden == 1:
        plt.ylabel('Concentración')
        plt.plot(tiempo, cinetica_orden_cero(tiempo, *parametros_optimizados[:2]), label='Ajuste de orden cero', color=color_linea)
        plt.legend().set_visible(False)
        plt.title('Orden 0')

    elif orden == 2:
        plt.ylabel('ln(Concentración)')
        plt.plot(tiempo, cinetica_orden_uno(tiempo, *parametros_optimizados[:2]), label='Ajuste de orden uno', color=color_linea)
        plt.legend().set_visible(False)
        plt.title('Orden 1')

    elif orden == 3:
        plt.ylabel('1 / Concentración')
        plt.plot(tiempo, cinetica_orden_dos(tiempo, *parametros_optimizados[:2]), label='Ajuste de orden dos', color=color_linea)
        plt.legend().set_visible(False)
        plt.title('Orden 2')

    k_cientifico = '{:.2e}'.format(parametros_optimizados[1])

    if orden == 1 or orden == 2:
        texto = f'k: {k_cientifico}\nOrdenada al origen: {parametros_optimizados[0]:.4f}\n$r^2$: {parametros_optimizados[2]:.4f}'
    else:
        texto = f'k: {k_cientifico}\nOrdenada al origen: {(1/parametros_optimizados[0]):.4f}\n$r^2$:{parametros_optimizados[2]:.4f}'

    plt.annotate(texto, xy=(0.3, 0.85), xycoords='axes fraction', fontsize=10, color='black')


# Ejemplo de uso
def analizar_cineticas(df, coef_extincion=1, longitud_camino=1):
    df = generar_dataframe_desde_Abs(df, coef_extincion, longitud_camino)

    # Inicializar listas para almacenar métricas de cada orden
    r2_values = []
    sse_values = []

    # Orden 0
    constante_ajuste_0, ordenada_al_origen_ajuste_0, r2_ajuste_0, sse_ajuste_0 = ajustar_cinetica_orden_cero(df['tiempo'], df['concentracion'])
    r2_values.append(r2_ajuste_0)
    sse_values.append(sse_ajuste_0)
    print("Orden 0:")
    print(f"Ordenada al origen de ajuste: {constante_ajuste_0}")
    print(f"Constante de ajuste: {ordenada_al_origen_ajuste_0}")
    print(f"Coeficiente de determinación (r²): {r2_ajuste_0}")

    # Orden 1
    constante_ajuste_1, ordenada_al_origen_ajuste_1, r2_ajuste_1, sse_ajuste_1 = ajustar_cinetica_orden_uno(df['tiempo'], np.log(df['concentracion']))
    r2_values.append(r2_ajuste_1)
    sse_values.append(sse_ajuste_1)
    print("\nOrden 1:")
    print(f"Ordenada al origen de ajuste: {constante_ajuste_1}")
    print(f"Constante de ajuste: {ordenada_al_origen_ajuste_1}")
    print(f"Coeficiente de determinación (r²): {r2_ajuste_1}")

    # Orden 2
    constante_ajuste_2, ordenada_al_origen_ajuste_2, r2_ajuste_2, sse_ajuste_2 = ajustar_cinetica_orden_dos(df['tiempo'], 1 / df['concentracion'])
    r2_values.append(r2_ajuste_2)
    sse_values.append(sse_ajuste_2)
    print("\nOrden 2:")
    print(f"Ordenada al origen de ajuste: {1/constante_ajuste_2}")
    print(f"Constante de ajuste: {ordenada_al_origen_ajuste_2}")
    print(f"Coeficiente de determinación (r²): {r2_ajuste_2}")

    # Gráficos
    plt.figure(figsize=(15, 5))

    graficar_cineticas(df['tiempo'], df['concentracion'], [constante_ajuste_0, ordenada_al_origen_ajuste_0, r2_ajuste_0], 1, color_puntos='blue', color_linea='blue')
    graficar_cineticas(df['tiempo'], np.log(df['concentracion']), [constante_ajuste_1, ordenada_al_origen_ajuste_1, r2_ajuste_1], 2, color_puntos='green', color_linea='green')
    graficar_cineticas(df['tiempo'], 1 / df['concentracion'], [constante_ajuste_2, ordenada_al_origen_ajuste_2, r2_ajuste_2], 3, color_puntos='red', color_linea='red')

    plt.tight_layout()
    plt.show()

    # Comparación de r2 y determinación del orden
    r2_values = [r2_ajuste_0, r2_ajuste_1, r2_ajuste_2]
    ordenes = ['Orden 0', 'Orden 1', 'Orden 2']
    mejor_orden = ordenes[np.argmax(r2_values)]

     # Determinar el mejor orden según r2 y SSE
    mejor_orden_r2 = np.argmax(r2_values)
    mejor_orden_sse = np.argmin(sse_values)

    print("\nResultados según r²:")
    print(f"Orden 0: {r2_values[0]}")
    print(f"Orden 1: {r2_values[1]}")
    print(f"Orden 2: {r2_values[2]}")

    print("\nResultados según SSE:")
    print(f"Orden 0: {sse_values[0]}")
    print(f"Orden 1: {sse_values[1]}")
    print(f"Orden 2: {sse_values[2]}")

    print(f"\nMejor orden según r²: Orden {mejor_orden_r2}")
    print(f"Mejor orden según SSE: Orden {mejor_orden_sse}")


def calcular_tiempo_vida_media(orden_reaccion, k, c0, t50):
    if orden_reaccion == 1:
        if k == 'no':
            k_calc = np.log(2) / t50
            t50_calc=t50
            c0_calc=c0
        if t50 == 'no':
            t50_calc = np.log(2) / k
            k_calc=k
            c0_calc=c0
    elif orden_reaccion == 2:
        if k == 'no':
            k_calc = 1 / (t50 * c0)
            t50_calc=t50
            c0_calc = c0
        if t50 == 'no':
            t50_calc = 1 / (k * c0)
            k_calc=k
            c0_calc = c0
        if c0 == 'no':
            c0_calc = 1 / (t50 * k)
            k_calc=k
            t50_calc=t50
    elif orden_reaccion == 0:
        if k == 'no':
            k_calc = c0 / (t50 * 2)
            t50_calc=t50
            c0_calc = c0
        if t50 == 'no':
            t50_calc = c0 / (2 * k)
            k_calc=k
            c0_calc = c0
        if c0 == 'no':
            c0 = t50 * 2 * k
            k_calc=k
            t50_calc=t50
    return k_calc, c0_calc, t50_calc

def calcular_tiempo_vida_util(orden_reaccion, k, c0, t90):
    if orden_reaccion == 1:
        if k == 'no':
            k_calc = np.log(0.9)*-1 / t90
            t90_calc=t90
            c0_calc=c0
        if t90 == 'no':
            t90_calc = np.log(0.9)*-1/ k
            k_calc=k
            c0_calc=c0
    elif orden_reaccion == 2:
        if k == 'no':
            k_calc = ((1/0.9)-1) / (t90 * c0)
            t90_calc=t90
            c0_calc = c0
        if t90 == 'no':
            t90_calc = ((1/0.9)-1) / (k * c0)
            k_calc=k
            c0_calc = c0
        if c0 == 'no':
            c0_calc = ((1/0.9)-1) / (t90 * k)
            k_calc=k
            t90_calc=t90
    elif orden_reaccion == 0:
        if k == 'no':
            k_calc = (c0*0.1) / (t90)
            t90_calc=t90
            c0_calc = c0
        if t90 == 'no':
            t90_calc = (c0*0.1) / (k)
            k_calc=k
            c0_calc = c0
        if c0 == 'no':
            c0 = (t90 * k)/0.1
            k_calc=k
            t90_calc=t90
    return k_calc, c0_calc, t90_calc


def ajustar_cinetica_orden_cero_multiple(tiempo, concentracion_triplicado):
    # Inicializar listas para almacenar los resultados de cada ajuste
    parametros_optimizados_lista = []
    r2_lista = []
    sse_lista = []

    # Iterar sobre cada set de datos de concentración
    for i in range(3):
        concentracion = concentracion_triplicado.iloc[:, i]

        # Ajuste para cada set de datos
        parametros_iniciales = [concentracion.iloc[0], 0.001]
        parametros_optimizados, matriz_covarianza = curve_fit(cinetica_orden_cero, tiempo, concentracion, p0=parametros_iniciales)

        concentracion_ajustada = cinetica_orden_cero(tiempo, *parametros_optimizados)
        r2 = r2_score(concentracion, concentracion_ajustada)
        sse = np.sum(np.square(concentracion - concentracion_ajustada))

        # Almacenar resultados en listas
        parametros_optimizados_lista.append(parametros_optimizados)
        r2_lista.append(r2)
        sse_lista.append(sse)

    # Calcular promedios y desviaciones estándar de parámetros
    parametros_optimizados_promedio = np.mean(parametros_optimizados_lista, axis=0)
    r2_promedio = np.mean(r2_lista)
    sse_promedio = np.mean(sse_lista)

    parametros_optimizados_std = np.std(parametros_optimizados_lista, axis=0)
    r2_std = np.std(r2_lista)
    sse_std = np.std(sse_lista)

    return parametros_optimizados_promedio, r2_promedio, sse_promedio, parametros_optimizados_std, r2_std, sse_std

def ajustar_cinetica_orden_uno_multiple(tiempo, concentracion_triplicado):
    # Inicializar listas para almacenar los resultados de cada ajuste
    parametros_optimizados_lista = []
    r2_lista = []
    sse_lista = []

    # Iterar sobre cada set de datos de concentración
    for i in range(3):
        concentracion_ln = np.log(concentracion_triplicado.iloc[:, i])

        # Ajuste para cada set de datos
        slope, intercept, r_value, p_value, std_err = linregress(tiempo, concentracion_ln)
        r2 = r_value**2
        sse = np.sum(np.square(concentracion_ln - (slope * tiempo + intercept)))

        # Almacenar resultados en listas
        parametros_optimizados_lista.append([intercept, -slope])
        r2_lista.append(r2)
        sse_lista.append(sse)

    # Calcular promedios y desviaciones estándar de parámetros
    parametros_optimizados_promedio = np.mean(parametros_optimizados_lista, axis=0)
    r2_promedio = np.mean(r2_lista)
    sse_promedio = np.mean(sse_lista)

    parametros_optimizados_std = np.std(parametros_optimizados_lista, axis=0)
    r2_std = np.std(r2_lista)
    sse_std = np.std(sse_lista)

    return parametros_optimizados_promedio, r2_promedio, sse_promedio, parametros_optimizados_std, r2_std, sse_std


def ajustar_cinetica_orden_dos_multiple(tiempo, inversa_concentracion_triplicado):
    inversa_concentracion_triplicado = 1 / inversa_concentracion_triplicado[['concentracion1', 'concentracion2', 'concentracion3']]
    # Inicializar listas para almacenar los resultados de cada ajuste
    c0_lista = []
    k_lista = []
    r2_lista = []
    sse_lista = []

    # Iterar sobre cada set de datos de inversa de concentración
    for i in range(3):
        inversa_concentracion = inversa_concentracion_triplicado.iloc[:, i]

        # Ajuste para cada set de datos
        parametros_iniciales = [1 / inversa_concentracion.iloc[0], 0.001]
        parametros_optimizados, matriz_covarianza = curve_fit(cinetica_orden_dos, tiempo, inversa_concentracion, p0=parametros_iniciales)

        concentracion_ajustada = cinetica_orden_dos(tiempo, *parametros_optimizados)
        r2 = r2_score(inversa_concentracion, concentracion_ajustada)
        sse = np.sum(np.square(inversa_concentracion - concentracion_ajustada))

        # Almacenar resultados en listas
        c0_lista.append(parametros_optimizados[0])
        k_lista.append(parametros_optimizados[1])
        r2_lista.append(r2)
        sse_lista.append(sse)

    # Calcular promedios y desviaciones estándar de parámetros
    c0_promedio = np.mean(c0_lista)
    k_promedio = np.mean(k_lista)
    r2_promedio = np.mean(r2_lista)
    sse_promedio = np.mean(sse_lista)

    c0_std = np.std(c0_lista)
    k_std = np.std(k_lista)
    r2_std = np.std(r2_lista)
    sse_std = np.std(sse_lista)

    parametros_optimizados_promedio = [1/c0_promedio, k_promedio]
    parametros_optimizados_std = [c0_std, k_std]

    return parametros_optimizados_promedio, r2_promedio, sse_promedio, parametros_optimizados_std, r2_std, sse_std

def graficar_cineticas_multiple(tiempo, datos, sd, parametros_optimizados_promedio, parametros_optimizados_std, r2, orden, color_puntos='blue', color_linea='red'):
    tiempo = np.array(tiempo)
    datos = np.array(datos)

    plt.subplot(1, 3, orden)
    sns.scatterplot(x=tiempo, y=datos, label='Datos experimentales', color=color_puntos)
    plt.xlabel('Tiempo')

    if orden == 1:
        plt.ylabel('Concentración')
        concentracion_ajustada = cinetica_orden_cero(tiempo, *parametros_optimizados_promedio[:2])
        plt.plot(tiempo, concentracion_ajustada, label='Ajuste de orden cero', color=color_linea)
        plt.legend().set_visible(False)
        #plt.fill_between(tiempo, df['AvSH13'] - std_SH13, df['AvSH13'] + std_SH13, color='red', alpha=0.2)
        plt.fill_between(tiempo,
                         concentracion_ajustada - sd,
                         concentracion_ajustada + sd,
                         color=color_linea, alpha=0.2)
        plt.title('Orden 0')

    elif orden == 2:
        plt.ylabel('ln(Concentración)')
        concentracion_ajustada = cinetica_orden_uno(tiempo, *parametros_optimizados_promedio[:2])
        plt.plot(tiempo, concentracion_ajustada, label='Ajuste de orden uno', color=color_linea)
        plt.legend().set_visible(False)
        plt.fill_between(tiempo,
                         concentracion_ajustada - sd,
                         concentracion_ajustada + sd,
                         color=color_linea, alpha=0.2)
        plt.title('Orden 1')

    elif orden == 3:
        plt.ylabel('1 / Concentración')
        c0 = parametros_optimizados_promedio[0]
        k = parametros_optimizados_promedio[1]
        inversa_concentracion_ajustada = cinetica_orden_dos(tiempo, 1/c0, k)
        
        plt.plot(tiempo, inversa_concentracion_ajustada, label='', color=color_linea)
        plt.legend().set_visible(False)
        plt.fill_between(tiempo,
                         inversa_concentracion_ajustada - sd,
                         inversa_concentracion_ajustada + sd,
                         color=color_linea, alpha=0.2)
        plt.title('Orden 2')

    k_cientifico = '{:.2e}'.format(parametros_optimizados_promedio[1])

    if orden == 1 or orden == 2:
        texto = f'k: {k_cientifico}\nOrdenada al origen: {parametros_optimizados_promedio[0]:.4f}\n$r^2$: {r2:.4f}'
    else:
        texto = f'k: {k_cientifico}\nOrdenada al origen: {(parametros_optimizados_promedio[0]):.4f}\n$r^2$:{r2:.4f}'

    plt.annotate(texto, xy=(0.3, 0.85), xycoords='axes fraction', fontsize=10, color='black')


def calcular_concentracion_multiple(absorbancia, coef_extincion=1, longitud_camino=1):
    concentracion = absorbancia / (longitud_camino * coef_extincion)
    return concentracion

def generar_dataframe_desde_Abs_multiple(df, coef_extincion=1, longitud_camino=1):
    if 'concentracion' not in df.columns:
        df['concentracion'] = calcular_concentracion_multiple(df['absorbancia'], coef_extincion, longitud_camino)
        df['ln_concentracion'] = np.log(df['concentracion'])
        df['inversa_concentracion'] = 1 / df['concentracion']
    print(df)
    return df

# Modificación de la función analizar_cineticas
def analizar_cineticas_multiple(df, coef_extincion=1, longitud_camino=1):
    #df = generar_dataframe_desde_Abs(df, coef_extincion, longitud_camino)

    # Inicializar listas para almacenar métricas de cada orden
    r2_values = []
    sse_values = []

    tiempo = df['tiempo']
    concentracion_triplicado = df[['concentracion1', 'concentracion2', 'concentracion3']]

    # Orden 0
    parametros_optimizados_promedio_0, r2_promedio_0, sse_promedio_0, parametros_optimizados_std_0, r2_std_0, sse_std_0 = ajustar_cinetica_orden_cero_multiple(tiempo, concentracion_triplicado)
    r2_values.append(r2_promedio_0)
    sse_values.append(sse_promedio_0)
    # Imprimir resultados
    print("\nResultados Promedio Orden 0:")
    print(f"Constante de ajuste: {parametros_optimizados_promedio_0[0]:.4f} ± {parametros_optimizados_std_0[0]:.4f}")
    print(f"Ordenada al origen: {parametros_optimizados_promedio_0[1]:.4f} ± {parametros_optimizados_std_0[1]:.4f}")
    print(f"Coeficiente de determinación (r²): {r2_promedio_0:.4f} ± {r2_std_0:.4f}")
    print(f"Suma de errores cuadráticos (SSE): {sse_promedio_0:.4f} ± {sse_std_0:.4f}")

    # Orden 1
    parametros_optimizados_promedio_1, r2_promedio_1, sse_promedio_1, parametros_optimizados_std_1, r2_std_1, sse_std_1 = ajustar_cinetica_orden_uno_multiple(tiempo, concentracion_triplicado)
    r2_values.append(r2_promedio_1)
    sse_values.append(sse_promedio_1)
    
    # Imprimir resultados
    print("\nResultados Promedio Orden 1:")
    print(f"Ordenada al origen: {parametros_optimizados_promedio_1[0]:.4f} ± {parametros_optimizados_std_1[0]:.4f}")
    print(f"Constante de ajuste (k): {parametros_optimizados_promedio_1[1]:.4f} ± {parametros_optimizados_std_1[1]:.4f}")
    print(f"Coeficiente de determinación (r²): {r2_promedio_1:.4f} ± {r2_std_1:.4f}")
    print(f"Suma de errores cuadráticos (SSE): {sse_promedio_1:.4f} ± {sse_std_1:.4f}")

    # Orden 2
    parametros_optimizados_promedio_2, r2_promedio_2, sse_promedio_2, parametros_optimizados_std_2, r2_std_2, sse_std_2 = ajustar_cinetica_orden_dos_multiple(tiempo, concentracion_triplicado)
    r2_values.append(r2_promedio_2)
    sse_values.append(sse_promedio_2)
    
    # Imprimir resultados
    print("\nResultados Promedio Orden 2:")
    print(f"Ordenada al origen: {parametros_optimizados_promedio_2[0]:.4f} ± {parametros_optimizados_std_2[0]:.4f}")
    print(f"Constante de ajuste (k): {parametros_optimizados_promedio_2[1]:.9f} ± {parametros_optimizados_std_2[1]:.4f}")
    print(f"Coeficiente de determinación (r²): {r2_promedio_2:.4f} ± {r2_std_2:.4f}")
    print(f"Suma de errores cuadráticos (SSE): {sse_promedio_2:.4f} ± {sse_std_2:.4f}")

    # Gráficos
    plt.figure(figsize=(15, 5))

    graficar_cineticas_multiple(df['tiempo'], df['concentracion_promedio'], df['std_concentracion'], parametros_optimizados_promedio_0, parametros_optimizados_std_0, r2_promedio_0, 1, color_puntos='blue', color_linea='blue')
    graficar_cineticas_multiple(df['tiempo'], np.log(df['concentracion_promedio']), df['std_log_concentracion'], parametros_optimizados_promedio_1, parametros_optimizados_std_1, r2_promedio_1, 2, color_puntos='green', color_linea='green')
    graficar_cineticas_multiple(df['tiempo'], 1 / df['concentracion_promedio'],  df['std_inversa_concentracion'], parametros_optimizados_promedio_2, parametros_optimizados_std_2, r2_promedio_2, 3, color_puntos='red', color_linea='red')

    plt.tight_layout()
    plt.show()

    # Comparación de r2 y determinación del orden
    r2_values = [r2_promedio_0, r2_promedio_1, r2_promedio_2]
    ordenes = ['Orden 0', 'Orden 1', 'Orden 2']
    mejor_orden_r2 = ordenes[np.argmax(r2_values)]
    mejor_orden_sse = ordenes[np.argmin(sse_values)]

    print("\nResultados según r²:")
    print(f"Orden 0: {r2_values[0]}")
    print(f"Orden 1: {r2_values[1]}")
    print(f"Orden 2: {r2_values[2]}")

    print("\nResultados según SSE:")
    print(f"Orden 0: {sse_values[0]}")
    print(f"Orden 1: {sse_values[1]}")
    print(f"Orden 2: {sse_values[2]}")

    print(f"\nMejor orden según r²: {mejor_orden_r2}")
    print(f"Mejor orden según SSE: {mejor_orden_sse}")
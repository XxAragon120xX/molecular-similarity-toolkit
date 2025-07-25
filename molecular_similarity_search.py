"""
Herramienta de Búsqueda por Similitud Molecular usando Coeficiente de Tanimoto
=============================================================================
Este script permite encontrar moléculas similares en una base de datos
basándose en la similitud estructural química utilizando fingerprints moleculares.
"""

# =============================================================================
# IMPORTACIONES
# =============================================================================
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
from PIL import Image

# RDKit para química computacional
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Draw, rdMolDescriptors, QED

# Configuración de warnings y logs
warnings.filterwarnings('ignore')
RDLogger.DisableLog('rdApp.*')

# =============================================================================
# CONFIGURACIÓN Y CONSTANTES
# =============================================================================
class Config:
    """Configuración centralizada del análisis."""
    
    # Rutas y archivos
    DATASET_URL = r'C:\Users\USUARIO\Desktop\Proyectos de Python\Autodiseñador\toxo.csv'
    OUTPUT_CSV_FOLDER = 'RT'
    OUTPUT_PNG_FOLDER = 'structures_png'
    
    # Parámetros de análisis
    SIMILARITY_THRESHOLD = 0.1
    NUM_MOLECULES_TO_VISUALIZE = 20
    FINGERPRINT_RADIUS = 2
    FINGERPRINT_BITS = 2048
    
    # Parámetros de visualización
    SUBIMG_SIZE = (250, 250)
    MOLS_PER_ROW = 3
    FIGURE_SIZE = (20, 20)
    
    # Nombres de archivos de salida
    @property
    def OUTPUT_CSV(self):
        return os.path.join(self.OUTPUT_CSV_FOLDER, 'filtered_compounds.csv')
    
    @property
    def OUTPUT_IMAGE(self):
        return os.path.join(self.OUTPUT_PNG_FOLDER, 'molecules.png')
    
    @property
    def OUTPUT_DISPLAYED_CSV(self):
        return os.path.join(self.OUTPUT_CSV_FOLDER, 'displayed_compounds.csv')

# =============================================================================
# CLASE PRINCIPAL PARA BÚSQUEDA DE SIMILITUD
# =============================================================================
class MolecularSimilaritySearch:
    """
    Clase principal para realizar búsquedas por similitud molecular.
    """
    
    def __init__(self, config=None):
        """
        Inicializa la herramienta de búsqueda por similitud.
        
        Parameters:
        -----------
        config : Config, optional
            Objeto de configuración. Si no se proporciona, usa Config por defecto.
        """
        self.config = config or Config()
        self.dataset = None
        self.query_molecule = None
        self.filtered_results = None
        
        # Crear carpetas de salida
        self._create_output_directories()
    
    def _create_output_directories(self):
        """Crea las carpetas de salida si no existen."""
        os.makedirs(self.config.OUTPUT_CSV_FOLDER, exist_ok=True)
        os.makedirs(self.config.OUTPUT_PNG_FOLDER, exist_ok=True)
    
    def load_dataset(self, dataset_path=None):
        """
        Carga y procesa el dataset molecular.
        
        Parameters:
        -----------
        dataset_path : str, optional
            Ruta al archivo del dataset. Si no se proporciona, usa la ruta por defecto.
        
        Returns:
        --------
        bool
            True si la carga fue exitosa, False en caso contrario.
        """
        try:
            path = dataset_path or self.config.DATASET_URL
            print(f"Cargando dataset desde: {path}")
            
            # Cargar datos
            self.dataset = pd.read_csv(path, delimiter=';')
            
            # Validar columnas requeridas
            required_columns = ['Smiles', 'Name', 'ChEMBL ID']
            missing_columns = [col for col in required_columns if col not in self.dataset.columns]
            
            if missing_columns:
                raise ValueError(f"Columnas faltantes: {missing_columns}")
            
            print(f"Dataset cargado: {len(self.dataset)} compuestos")
            
            # Procesar SMILES
            self._process_smiles()
            
            print(f"Compuestos válidos después del procesamiento: {len(self.dataset)}")
            return True
            
        except Exception as e:
            print(f"Error al cargar dataset: {e}")
            return False
    
    def _process_smiles(self):
        """Procesa las estructuras SMILES y crea objetos moleculares."""
        print("Procesando estructuras SMILES...")
        
        # Limpiar y validar SMILES
        self.dataset['Smiles'] = self.dataset['Smiles'].apply(
            lambda x: str(x) if pd.notnull(x) else None
        )
        
        # Convertir SMILES a objetos moleculares
        self.dataset['Structure'] = self.dataset['Smiles'].apply(
            lambda x: Chem.MolFromSmiles(x) if x else None
        )
        
        # Eliminar moléculas inválidas
        initial_count = len(self.dataset)
        self.dataset = self.dataset.dropna(subset=['Structure'])
        
        removed_count = initial_count - len(self.dataset)
        if removed_count > 0:
            print(f"Se removieron {removed_count} estructuras inválidas")
    
    def set_query_molecule(self, smiles_string):
        """
        Establece la molécula de consulta a partir de un SMILES.
        
        Parameters:
        -----------
        smiles_string : str
            String SMILES de la molécula de consulta.
        
        Returns:
        --------
        bool
            True si la molécula es válida, False en caso contrario.
        """
        try:
            self.query_molecule = Chem.MolFromSmiles(smiles_string)
            
            if self.query_molecule is None:
                print("Error: No se pudo crear la molécula a partir del SMILES proporcionado")
                return False
            
            print(f"Molécula de consulta establecida: {smiles_string}")
            return True
            
        except Exception as e:
            print(f"Error al procesar SMILES de consulta: {e}")
            return False
    
    def calculate_similarities(self):
        """
        Calcula las similitudes de Tanimoto entre la molécula de consulta y el dataset.
        
        Returns:
        --------
        list
            Lista de valores de similitud de Tanimoto.
        """
        if self.query_molecule is None:
            raise ValueError("No se ha establecido una molécula de consulta")
        
        if self.dataset is None:
            raise ValueError("No se ha cargado un dataset")
        
        print("Calculando similitudes de Tanimoto...")
        
        # Generar fingerprint de la molécula de consulta
        query_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
            self.query_molecule, 
            self.config.FINGERPRINT_RADIUS, 
            nBits=self.config.FINGERPRINT_BITS
        )
        
        # Generar fingerprints del dataset
        dataset_fps = [
            rdMolDescriptors.GetMorganFingerprintAsBitVect(
                mol, 
                self.config.FINGERPRINT_RADIUS, 
                nBits=self.config.FINGERPRINT_BITS
            ) 
            for mol in self.dataset['Structure']
        ]
        
        # Calcular similitudes
        similarities = [
            DataStructs.FingerprintSimilarity(query_fp, fp) 
            for fp in dataset_fps
        ]
        
        print(f"Similitudes calculadas para {len(similarities)} compuestos")
        return similarities
    
    def filter_and_rank_results(self, similarities, threshold=None):
        """
        Filtra y ordena los resultados basándose en la similitud.
        
        Parameters:
        -----------
        similarities : list
            Lista de valores de similitud.
        threshold : float, optional
            Umbral mínimo de similitud. Si no se proporciona, usa el valor por defecto.
        
        Returns:
        --------
        pd.DataFrame
            DataFrame filtrado y ordenado.
        """
        threshold = threshold or self.config.SIMILARITY_THRESHOLD
        
        # Añadir similitudes al dataset
        dataset_copy = self.dataset.copy()
        dataset_copy['tanimoto_similarity'] = similarities
        
        # Filtrar por umbral
        filtered = dataset_copy[dataset_copy['tanimoto_similarity'] >= threshold]
        
        # Ordenar por similitud descendente
        self.filtered_results = filtered.sort_values(
            by='tanimoto_similarity', 
            ascending=False
        ).reset_index(drop=True)
        
        print(f"Compuestos que superan el umbral ({threshold}): {len(self.filtered_results)}")
        
        return self.filtered_results
    
    def save_results(self, save_path=None):
        """
        Guarda los resultados filtrados en un archivo CSV.
        
        Parameters:
        -----------
        save_path : str, optional
            Ruta donde guardar el archivo. Si no se proporciona, usa la ruta por defecto.
        """
        if self.filtered_results is None:
            print("No hay resultados para guardar")
            return
        
        save_path = save_path or self.config.OUTPUT_CSV
        self.filtered_results.to_csv(save_path, index=False)
        print(f"Resultados guardados en: {save_path}")
    
    def visualize_top_molecules(self, num_molecules=None, save_image=True):
        """
        Visualiza las moléculas más similares.
        
        Parameters:
        -----------
        num_molecules : int, optional
            Número de moléculas a visualizar. Si no se proporciona, usa el valor por defecto.
        save_image : bool
            Si True, guarda la imagen generada.
        
        Returns:
        --------
        PIL.Image
            Imagen con la grilla de moléculas.
        """
        if self.filtered_results is None or len(self.filtered_results) == 0:
            print("No hay resultados para visualizar")
            return None
        
        num_molecules = num_molecules or self.config.NUM_MOLECULES_TO_VISUALIZE
        
        # Preparar datos para visualización
        molecules_to_show = []
        legends = []
        displayed_data = []
        
        for i in range(min(num_molecules, len(self.filtered_results))):
            row = self.filtered_results.iloc[i]
            mol = row['Structure']
            
            if mol is not None:
                molecules_to_show.append(mol)
                
                # Crear leyenda informativa
                legend = (
                    f"{row['Name']}\n"
                    f"Tanimoto: {row['tanimoto_similarity']:.3f}\n"
                    f"ChEMBL ID: {row['ChEMBL ID']}"
                )
                legends.append(legend)
                displayed_data.append(row)
        
        if not molecules_to_show:
            print("No hay moléculas válidas para visualizar")
            return None
        
        print(f"Visualizando {len(molecules_to_show)} moléculas")
        
        # Generar imagen
        img = Draw.MolsToGridImage(
            molecules_to_show,
            legends=legends,
            subImgSize=self.config.SUBIMG_SIZE,
            molsPerRow=self.config.MOLS_PER_ROW
        )
        
        # Guardar información de moléculas mostradas
        displayed_df = pd.DataFrame(displayed_data)
        displayed_df.to_csv(self.config.OUTPUT_DISPLAYED_CSV, index=False)
        print(f"Información de moléculas mostradas guardada en: {self.config.OUTPUT_DISPLAYED_CSV}")
        
        # Guardar imagen si se solicita
        if save_image:
            img.save(self.config.OUTPUT_IMAGE)
            print(f"Imagen guardada en: {self.config.OUTPUT_IMAGE}")
        
        return img
    
    def display_results(self, img):
        """
        Muestra los resultados de la visualización.
        
        Parameters:
        -----------
        img : PIL.Image
            Imagen a mostrar.
        """
        if img is None:
            return
        
        plt.figure(figsize=self.config.FIGURE_SIZE)
        plt.imshow(img)
        plt.axis('off')
        plt.title('Moléculas Más Similares', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.show()
    
    def get_summary_statistics(self):
        """
        Proporciona estadísticas resumen de los resultados.
        
        Returns:
        --------
        dict
            Diccionario con estadísticas resumen.
        """
        if self.filtered_results is None:
            return {}
        
        similarities = self.filtered_results['tanimoto_similarity']
        
        stats = {
            'total_compounds': len(self.filtered_results),
            'max_similarity': similarities.max(),
            'min_similarity': similarities.min(),
            'mean_similarity': similarities.mean(),
            'median_similarity': similarities.median(),
            'std_similarity': similarities.std()
        }
        
        return stats

# =============================================================================
# FUNCIÓN PRINCIPAL INTERACTIVA
# =============================================================================
def main():
    """
    Función principal que ejecuta el análisis de similitud molecular de forma interactiva.
    """
    print("=" * 60)
    print("HERRAMIENTA DE BÚSQUEDA POR SIMILITUD MOLECULAR")
    print("=" * 60)
    
    # Inicializar herramienta
    similarity_search = MolecularSimilaritySearch()
    
    # Cargar dataset
    if not similarity_search.load_dataset():
        print("Error al cargar el dataset. Terminando...")
        return
    
    # Solicitar molécula de consulta
    print("\n" + "-" * 40)
    query_smiles = input("Ingrese el SMILES de la molécula de consulta: ").strip()
    
    if not similarity_search.set_query_molecule(query_smiles):
        print("Error con la molécula de consulta. Terminando...")
        return
    
    # Calcular similitudes
    print("\n" + "-" * 40)
    similarities = similarity_search.calculate_similarities()
    
    # Filtrar y ordenar resultados
    filtered_results = similarity_search.filter_and_rank_results(similarities)
    
    if len(filtered_results) == 0:
        print(f"No se encontraron compuestos con similitud >= {similarity_search.config.SIMILARITY_THRESHOLD}")
        return
    
    # Mostrar estadísticas
    print("\n" + "-" * 40)
    print("ESTADÍSTICAS DE RESULTADOS:")
    stats = similarity_search.get_summary_statistics()
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"{key}: {value:.4f}")
        else:
            print(f"{key}: {value}")
    
    # Guardar resultados
    print("\n" + "-" * 40)
    similarity_search.save_results()
    
    # Visualizar resultados
    print("Generando visualización...")
    img = similarity_search.visualize_top_molecules()
    
    if img:
        similarity_search.display_results(img)
    
    print("\n" + "=" * 60)
    print("ANÁLISIS COMPLETADO")
    print("=" * 60)
    print(f"Archivos generados:")
    print(f"- {similarity_search.config.OUTPUT_CSV}")
    print(f"- {similarity_search.config.OUTPUT_DISPLAYED_CSV}")
    print(f"- {similarity_search.config.OUTPUT_IMAGE}")

# =============================================================================
# FUNCIONES DE UTILIDAD ADICIONALES
# =============================================================================
def batch_similarity_search(query_smiles_list, dataset_path=None, threshold=0.1):
    """
    Realiza búsquedas de similitud para múltiples moléculas de consulta.
    
    Parameters:
    -----------
    query_smiles_list : list
        Lista de SMILES para buscar.
    dataset_path : str, optional
        Ruta al dataset.
    threshold : float
        Umbral de similitud.
    
    Returns:
    --------
    dict
        Diccionario con resultados para cada consulta.
    """
    results = {}
    similarity_search = MolecularSimilaritySearch()
    
    # Cargar dataset una sola vez
    if not similarity_search.load_dataset(dataset_path):
        return results
    
    for i, query_smiles in enumerate(query_smiles_list):
        print(f"\nProcesando consulta {i+1}/{len(query_smiles_list)}: {query_smiles}")
        
        if similarity_search.set_query_molecule(query_smiles):
            similarities = similarity_search.calculate_similarities()
            filtered_results = similarity_search.filter_and_rank_results(similarities, threshold)
            results[query_smiles] = filtered_results
        else:
            results[query_smiles] = None
    
    return results

# =============================================================================
# EJECUCIÓN
# =============================================================================
if __name__ == "__main__":
    main()

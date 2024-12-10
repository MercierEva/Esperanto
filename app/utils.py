import os
import time
import gzip
import shutil

def get_file_info(file_path):
    """Retourne les informations basiques d'un fichier: nom, extension, taille, date de modification."""
    try:
        basename = os.path.basename(file_path)
        name, ext = os.path.splitext(basename)
        size = os.path.getsize(file_path)
        modification_time = os.path.getmtime(file_path)
        return name, ext, size, modification_time
    except FileNotFoundError:
        raise FileNotFoundError(f"Le fichier {file_path} n'a pas été trouvé.")
    except Exception as e:
        raise Exception(f"Erreur lors de la récupération des informations du fichier: {str(e)}")

def format_size(size):
    """Formate la taille du fichier en Ko, Mo, etc."""
    if size < 1024:
        return f"{size} B"
    elif size < 1024**2:
        return f"{size / 1024:.2f} KB"
    elif size < 1024**3:
        return f"{size / 1024**2:.2f} MB"
    else:
        return f"{size / 1024**3:.2f} GB"

def get_file_modification_time(file_path):
    """Retourne la date de dernière modification d'un fichier sous forme lisible."""
    try:
        mod_time = os.path.getmtime(file_path)
        return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(mod_time))
    except FileNotFoundError:
        raise FileNotFoundError(f"Le fichier {file_path} n'a pas été trouvé.")
    except Exception as e:
        raise Exception(f"Erreur lors de la récupération de la date de modification: {str(e)}")

def read_gzip_file(file_path):
    """Lit un fichier gzip et retourne son contenu sous forme de liste de lignes."""
    try:
        with gzip.open(file_path, 'rt') as file:
            return file.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"Le fichier {file_path} n'a pas été trouvé.")
    except Exception as e:
        raise Exception(f"Erreur lors de la lecture du fichier gzip: {str(e)}")

def delete_files_in_directory(directory_path):
    """Supprime tous les fichiers dans un répertoire donné."""
    try:
        if os.path.exists(directory_path):
            for file_name in os.listdir(directory_path):
                file_path = os.path.join(directory_path, file_name)
                if os.path.isfile(file_path):
                    os.remove(file_path)
    except Exception as e:
        raise Exception(f"Erreur lors de la suppression des fichiers: {str(e)}")

def copy_file(src, dst):
    """Copie un fichier d'un emplacement à un autre."""
    try:
        shutil.copy(src, dst)
    except FileNotFoundError:
        raise FileNotFoundError(f"Le fichier source {src} n'a pas été trouvé.")
    except Exception as e:
        raise Exception(f"Erreur lors de la copie du fichier: {str(e)}")

def get_files_in_directory(directory_path, extension=None):
    """Retourne une liste de fichiers dans un répertoire, avec une option de filtrage par extension."""
    try:
        files = []
        for file_name in os.listdir(directory_path):
            file_path = os.path.join(directory_path, file_name)
            if os.path.isfile(file_path):
                if extension is None or file_name.endswith(extension):
                    files.append(file_path)
        return files
    except Exception as e:
        raise Exception(f"Erreur lors de la récupération des fichiers: {str(e)}")

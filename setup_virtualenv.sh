#!/bin/bash

# Fonction pour installer les dépendances système
install_dependencies() {
    # Vérifier si le système est Ubuntu ou Fedora
    if [ -f /etc/lsb-release ]; then
        # Ubuntu/Debian
        echo "Installation des dépendances pour Ubuntu/Debian..."
        sudo apt update
        sudo apt install -y build-essential dpkg-dev python3-dev freeglut3-dev \
            libgl1-mesa-dev libglu1-mesa-dev libgstreamer-plugins-base1.0-dev \
            libgtk-3-dev libjpeg-dev libnotify-dev libpng-dev libsdl2-dev \
            libsm-dev libtiff-dev libwebkit2gtk-4.0-dev libxtst-dev
    elif [ -f /etc/fedora-release ]; then
        # Fedora
        echo "Installation des dépendances pour Fedora..."
        sudo dnf update -y
        sudo dnf groupinstall -y "Development Tools"
        sudo dnf install -y python3-devel freeglut-devel mesa-libGL-devel \
            mesa-libGLU-devel gstreamer1-plugins-base-devel gtk3-devel \
            libjpeg-devel libnotify-devel libpng-devel SDL2-devel libSM-devel \
            libtiff-devel webkit2gtk3-devel libXtst-devel
    else
        echo "Ce script est conçu pour Ubuntu/Debian et Fedora uniquement."
        exit 1
    fi
}

# Fonction pour configurer un environnement virtuel avec virtualenv
setup_virtualenv() {
    # Installer virtualenv si nécessaire
    echo "Installation de virtualenv..."
    python3 -m pip install --upgrade pip
    python3 -m pip install virtualenv

    # Créer un environnement virtuel
    echo "Création de l'environnement virtuel 'env_wxpython'..."
    python3 -m venv env_wxpython

    # Activer l'environnement virtuel
    source env_wxpython/bin/activate
}

# Fonction pour installer wxPython dans l'environnement virtuel
install_wxpython() {
    # Installer wxPython
    echo "Installation de wxPython..."
    pip install wxPython pandas snakemake PyYAML
}

# Exécution des fonctions
install_dependencies
setup_virtualenv
install_wxpython

echo "Configuration terminée. L'environnement virtuel 'env_wxpython' a été créé et wxPython a été installé."


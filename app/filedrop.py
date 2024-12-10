import wx
import os

class FileDrop(wx.FileDropTarget):
    """
    Classe pour gérer l'importation de fichiers via un drag-and-drop.
    """
    
    def __init__(self, listctrl):
        """
        Constructeur de la classe FileDrop.
        :param listctrl: Liste de contrôle où les fichiers seront ajoutés après avoir été déposés.
        """
        super().__init__()
        self.listctrl = listctrl  # Référence au contrôle de liste (ListCtrl) où les fichiers seront ajoutés

    def OnDropFiles(self, x, y, filenames):
        """
        Méthode appelée lorsque des fichiers sont déposés.
        :param x, y: Coordonnées de la position du curseur lorsque les fichiers sont déposés.
        :param filenames: Liste des chemins des fichiers déposés.
        :return: True si l'opération de dépôt a réussi.
        """
        for file in filenames:
            self.listctrl.updateText(file)  # Met à jour le contrôle de liste avec les fichiers déposés
        return True

import os
import time
import wx

class MyListCtrl(wx.ListCtrl):
    """
    Classe personnalisée héritée de wx.ListCtrl pour afficher des informations
    sur les fichiers dans une liste avec des colonnes spécifiques.
    """

    def __init__(self, parent, id):
        """Initialise le contrôle de la liste avec des colonnes spécifiques."""
        super().__init__(parent, id, style=wx.LC_REPORT)
        self.c = 0

        # Définition des colonnes
        self.InsertColumn(0, 'Name')
        self.InsertColumn(1, 'Ext')
        self.InsertColumn(2, 'Size', wx.LIST_FORMAT_RIGHT)
        self.InsertColumn(3, 'Modified')
        self.InsertColumn(4, 'Path')

        # Largeur des colonnes
        self.SetColumnWidth(0, 130)
        self.SetColumnWidth(1, 70)
        self.SetColumnWidth(2, 100)
        self.SetColumnWidth(3, 130)
        self.SetColumnWidth(4, 70)

    def update_text(self, file_path):
        """Met à jour le contrôle de la liste avec les informations d'un fichier."""
        basename = os.path.basename(file_path)
        name, ext = os.path.splitext(basename)
        name_cut, extcomp = os.path.splitext(name)
        ex = extcomp + ext
        size = os.path.getsize(file_path)
        modification_time = os.path.getmtime(file_path)

        # Insertion de l'élément dans la liste
        self.InsertItem(self.c, name_cut)
        self.SetItem(self.c, 1, ex)
        self.SetItem(self.c, 2, f"{size} B")
        self.SetItem(self.c, 3, time.strftime('%Y-%m-%d %H:%M', time.localtime(modification_time)))
        self.SetItem(self.c, 4, file_path)

        # Couleur de fond alternative pour les lignes
        if self.c % 2 == 0:
            self.SetItemBackgroundColour(self.c, '#e6f1f5')

        self.c += 1

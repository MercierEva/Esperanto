import os
import sys
import time
import glob
import gzip
import queue
import threading
import subprocess
import wx
from wx import lib
from wx.lib.pubsub import pub
from utils import get_file_info, read_gzip_file, format_size, get_file_modification_time, delete_files_in_directory
from listctrl import MyListCtrl
from filedrop import FileDropTarget
from dialogs import ExceptionDialog, ConfirmationDialog


class MyFrame(wx.Frame):
    """
    Classe principale qui définit l'interface et les fonctionnalités principales de l'application.
    """

    def __init__(self, parent, title):
        super().__init__(parent, title=title, size=(1300, 700))
        self.workspace = os.getcwd()
        self.init_ui()

    def init_ui(self):
        """Initialise l'interface utilisateur."""
        self.setup_icon()
        self.listbook = self.create_listbook()
        self.setup_statusbar()
        self.bind_events()
        self.setup_filedrop()

    def setup_icon(self):
        """Configure l'icône de la fenêtre."""
        styleicon_path = os.path.join(self.workspace, "resources/images/clipart1150338.png")
        styleicon = wx.Image(styleicon_path)
        styleicon.Rescale(20, 20)
        self.SetIcon(wx.Icon(styleicon.ConvertToBitmap()))

    def create_listbook(self):
        """Crée le composant Listbook et configure ses options."""
        listbook = Listbook(self)
        listbook.Children[0].SetBackgroundColour('#e6e6ff')
        listbook.SetBackgroundColour('#93c2cf')
        return listbook

    def setup_statusbar(self):
        """Configure la barre de statut avec une barre de progression."""
        self.statusbar = self.CreateStatusBar(3)
        self.text1 = wx.StaticText(self.statusbar, wx.ID_ANY, style=wx.ALIGN_LEFT)
        self.progress_bar = wx.Gauge(self.statusbar, wx.ID_ANY, style=wx.GA_HORIZONTAL | wx.GA_SMOOTH | wx.ALIGN_RIGHT)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.text1, 1, wx.EXPAND | wx.ALL, 5)
        sizer.Add(self.progress_bar, 1, wx.EXPAND | wx.LEFT | wx.BOTTOM, 5)
        self.statusbar.SetSizer(sizer)
        self.progress_bar.SetRange(100)
        self.progress_bar.SetValue(0)

    def bind_events(self):
        """Lie les événements des boutons aux méthodes correspondantes."""
        self.listbook.buttonexe1.Bind(wx.EVT_BUTTON, self.generate_histogram)
        self.listbook.buttonRun.Bind(wx.EVT_BUTTON, self.run_snakemake)
        self.listbook.killbutton.Bind(wx.EVT_BUTTON, self.kill_process)

    def setup_filedrop(self):
        """Configure le gestionnaire de dépôt de fichiers."""
        self.file_drop_target = FileDropTarget(self.listbook.listarr0)
        self.listbook.listarr0.SetDropTarget(self.file_drop_target)

    def generate_histogram(self, event):
        """Génère un histogramme en fonction des fichiers ajoutés."""
        self.text1.SetLabelText("Started")
        self.progress_bar.Show()
        self.progress_bar.Pulse()
        self.listbook.SetSelection(1)

        dico_data_input = self.collect_file_data()
        self.process_file_data(dico_data_input)

        self.text1.SetLabelText("Finished")
        self.progress_bar.Hide()

    def collect_file_data(self):
        """Collecte les informations des fichiers à partir du Listbook."""
        dico_data_input = {}
        for i in range(self.listbook.listarr0.GetItemCount()):
            name_file = self.listbook.listarr0.GetItemText(item=i, col=0)
            pathfile = self.listbook.listarr0.GetItemText(item=i, col=4)
            dico_data_input[name_file] = pathfile
        return dico_data_input

    def process_file_data(self, dico_data_input):
        """Lit les fichiers et génère des données d'histogrammes."""
        lens = []
        for name, path in dico_data_input.items():
            self.listbook.TextCtrlRunHisto.AppendText(f"{name}\r\n")
            data = read_gzip_file(path)
            lens.extend(len(data[l].strip()) for l in range(1, len(data) // 4, 4))
        self.generate_histogram_output(sorted(lens))

    def generate_histogram_output(self, lens_s):
        """Génère l'affichage de l'histogramme."""
        step = 100
        max_range = int(self.listbook.interval)
        for k in range(0, max_range, step):
            self.listbook.TextCtrlRunHisto.AppendText(f"{k}-{k + step}  ")
            count = len([x for x in lens_s if k <= x < k + step])
            percentage = count * 100 // len(lens_s)
            self.listbook.TextCtrlRunHisto.AppendText(percentage * 'X')
            self.listbook.TextCtrlRunHisto.AppendText(f" {percentage}%\n")

    def kill_process(self, event):
        """Interrompt le processus en cours."""
        if hasattr(self, 'process_sk') and self.process_sk:
            self.process_sk.terminate()
            self.process_sk.returncode = -1
            self.text1.SetLabelText("Process killed!")

    def run_snakemake(self, event):
        """Exécute Snakemake."""
        try:
            self.text1.SetLabelText("Started")
            self.progress_bar.Show()

            self.clean_snakemake_temp_files()
            self.listbook.fill_config_snakemake()

            self.process_sk = subprocess.Popen(
                f"snakemake -s {self.workspace}/workflow/SnakeFile --cores {self.listbook.thread} --use-conda",
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

            self.monitor_process_output()
            self.text1.SetLabelText("Finished")
            self.progress_bar.Hide()

        except ValueError:
            dlg = ExceptionDialog("Please fill in the fields correctly!")
            dlg.ShowModal()

    def clean_snakemake_temp_files(self):
        """Nettoie les fichiers temporaires de Snakemake."""
        delete_files_in_directory('.snakemake/locks')
        delete_files_in_directory('.snakemake/incomplete')

    def monitor_process_output(self):
        """Surveille la sortie du processus Snakemake."""
        pa_q = queue.Queue()
        pa_t = threading.Thread(target=self.read_output, args=(self.process_sk.stdout, pa_q))
        pa_t.daemon = True
        pa_t.start()

        while True:
            self.progress_bar.Pulse()
            wx.GetApp().Yield()
            self.process_sk.poll()
            try:
                line = pa_q.get(False)
                self.listbook.TextCtrlRunSK.AppendText(f"\n{line}")
            except queue.Empty:
                pass
            if self.process_sk.returncode is not None:
                break

    def read_output(self, pipe, q):
        """Lit les sorties du processus et les transmet à une file d'attente."""
        for line in iter(pipe.readline, b''):
            q.put(line.decode())


def main():
    app = wx.App(False)  # False pour ne pas initialiser un redirector d'affichage (console)
    ex = MyFrame(None, 'Esperanto')  # Créer une instance de la classe MyFrame
    ex.Show()  # Afficher la fenêtre principale
    app.MainLoop()  # Démarrer la boucle principale de l'application

if __name__ == '__main__':
    main()  # Appeler la fonction main si ce fichier est exécuté directement

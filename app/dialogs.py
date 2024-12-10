import wx
import wx.lib.dialogs as GMD  # Importation de la bibliothèque de dialogues génériques

class ExceptionDialog(GMD.GenericMessageDialog):
    """
    Classe pour gérer les boîtes de dialogue d'exception.
    Hérite de wx.lib.dialogs.GenericMessageDialog pour afficher une erreur sous forme de message générique.
    """

    def __init__(self, msg):
        """
        Constructeur pour initialiser le message et la boîte de dialogue.
        :param msg: Le message d'erreur à afficher.
        """
        super().__init__(None, msg, "Exception!", wx.OK | wx.ICON_ERROR)

class InputDialog(wx.Dialog):
    """
    Classe pour gérer une boîte de dialogue d'entrée personnalisée.
    Permet à l'utilisateur de saisir des informations.
    """

    def __init__(self, parent, title, message):
        """
        Constructeur de la boîte de dialogue d'entrée.
        :param parent: Fenêtre parente.
        :param title: Titre de la boîte de dialogue.
        :param message: Message à afficher dans la boîte de dialogue.
        """
        super().__init__(parent, title=title, size=(300, 200))
        
        # Crée un panneau pour la boîte de dialogue
        panel = wx.Panel(self)
        self.sizer = wx.BoxSizer(wx.VERTICAL)

        # Ajout du message
        self.msg = wx.StaticText(panel, label=message)
        self.sizer.Add(self.msg, 0, wx.ALL | wx.EXPAND, 10)

        # Zone de texte pour l'entrée utilisateur
        self.input_text = wx.TextCtrl(panel)
        self.sizer.Add(self.input_text, 0, wx.ALL | wx.EXPAND, 10)

        # Boutons OK et Cancel
        button_sizer = self.CreateButtonSizer(wx.OK | wx.CANCEL)
        self.sizer.Add(button_sizer, 0, wx.ALL | wx.EXPAND, 5)

        panel.SetSizerAndFit(self.sizer)

    def get_input(self):
        """
        Retourne le texte saisi par l'utilisateur.
        """
        return self.input_text.GetValue()

class ConfirmationDialog(wx.MessageDialog):
    """
    Classe pour gérer les boîtes de dialogue de confirmation.
    """

    def __init__(self, parent, message, caption="Confirmation"):
        """
        Constructeur pour afficher un message de confirmation.
        :param parent: Fenêtre parente.
        :param message: Message à afficher.
        :param caption: Titre de la boîte de dialogue.
        """
        super().__init__(parent, message, caption, wx.YES_NO | wx.ICON_QUESTION)

    def get_response(self):
        """
        Retourne True si l'utilisateur a cliqué sur 'Oui', False sinon.
        """
        return self.GetReturnCode() == wx.YES

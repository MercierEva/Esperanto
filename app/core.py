#!/bin/env python
#-*- coding: utf-8 -*-

import wx
from wx import html
from wx.html2 import *
import sys
import traceback
import pathlib
import os
import subprocess
import yaml
import wx.lib.agw.genericmessagedialog as GMD
import time
import re
import gzip
import glob
import threading
import queue
from wx import lib
from pubsub import pub
from .utils import Listbook
from .listctrl import MyListCtrl
from .filedrop import FileDrop


path_to_home = os.path.expanduser('~')
path_to_local = os.path.join(path_to_home + '/.local')
path_to_gtk = os.path.join(path_to_local +'/share')

access_rights = 0o755

try:
    os.mkdir(path_to_local, access_rights)
    os.mkdir(path_to_gtk, access_rights)
except FileExistsError:
    pass



class ExceptionDialog(GMD.GenericMessageDialog):
    """"""
    #----------------------------------------------------------------------
    def __init__(self, msg):
        """Constructor"""
        GMD.GenericMessageDialog.__init__(self, None, msg, "Exception!",
                                          wx.OK|wx.ICON_ERROR)

def MyExceptionHook(etype, value, trace):
    """
    Handler for all unhandled exceptions.
    :param `etype`: the exception type (`SyntaxError`, `ZeroDivisionError`, etc...);
    :type `etype`: `Exception`
    :param string `value`: the exception error message;
    :param string `trace`: the traceback header, if any (otherwise, it prints the
    standard Python header: ``Traceback (most recent call last)``.
    """
    frame = wx.GetApp().GetTopWindow()
    tmp = traceback.format_exception(etype, value, trace)
    exception = "".join(tmp)
    
    dlg = ExceptionDialog(exception)
    dlg.ShowModal()
    if dlg.ShowModal == wx.ID_OK or wx.ID_CANCEL:
        dlg.Destroy()    



class MyFrame(wx.Frame):

    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title,  size=(1300,700))
        sys.excepthook = MyExceptionHook
        self.workspace = os.getcwd()
        styleicon=wx.Image(self.workspace + "/images/clipart1150338.png")
        styleicon.Rescale(20, 20)
        icon = styleicon.ConvertToBitmap()

        self.SetIcon(wx.Icon(icon)) 
        self.listbook = Listbook(self)
        self.listbook.Children[0].SetBackgroundColour('#e6e6ff') 

        self.listbook.SetBackgroundColour('#93c2cf')
        self.listbook.buttonexe1.Bind(wx.EVT_BUTTON, self.generateHisto)
        self.listbook.buttonRun.Bind(wx.EVT_BUTTON, self.Running_Snakemake)
        self.listbook.killbutton.Bind(wx.EVT_BUTTON, self.kill_process)

        self.statusbar = self.CreateStatusBar(3)
        self.text1 = wx.StaticText(self.statusbar, wx.ID_ANY, style=wx.ALIGN_LEFT)
        self.progress_bar = wx.Gauge(self.statusbar, wx.ID_ANY, style=wx.GA_HORIZONTAL | wx.GA_SMOOTH | wx.ALIGN_RIGHT)
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer.Add(self.text1, 1, wx.EXPAND | wx.ALL, 5)
        self.sizer.Add(self.progress_bar, 1, wx.EXPAND | wx.LEFT | wx.BOTTOM, 5)
        self.statusbar.SetSizer(self.sizer)
        self.progress_bar.SetRange(100)
        self.progress_bar.SetValue(0)

    
    def generateHisto(self, event):
        self.text1.SetLabelText("Started")     
        self.progress_bar.Show()  
        self.progress_bar.Pulse()
        self.listbook.SetSelection(1) 
 
        dico_data_input = {}
        i=0
        while i < self.listbook.listarr0.GetItemCount():
            name_file = self.listbook.listarr0.GetItemText(item=i, col=0)
            pathfile = self.listbook.listarr0.GetItemText(item=i, col=4)
            dico_data_input[str(name_file)] = str(pathfile)  
            i += 1
 
        lens=[]

        for k,v in dico_data_input.items() :
            self.listbook.TextCtrlRunHisto.AppendText(k + "\r\n")
            with gzip.open(v, 'r') as fq:
                data = fq.readlines()   
                for l in range(1, len(data)//4, 4):
                    lens.append(len(data[l].strip()))
            lens_s = sorted(lens) 
            step = 100
            max_range=int(self.listbook.interval)
            for k in range(0, max_range, step):
                self.listbook.TextCtrlRunHisto.AppendText("%d-%d  " % (k, k+ step))
                j = filter(lambda x: x >= k and x < k+step, lens_s)
                p = len(list(j))*100//len(lens_s)
                self.listbook.TextCtrlRunHisto.AppendText( p * 'X')
                self.listbook.TextCtrlRunHisto.AppendText(' % \n')

            self.text1.SetLabelText("Finished")
            self.progress_bar.Hide()
            

    def kill_process(self, event):
        self.process_sk.terminate()
        self.process_sk.returncode = -1 
        self.text1.SetLabelText("Process killed !")
    
    def read_output(self, pipe, q):
        while True:
            l = pipe.readline() 
            if l :
                newoutput = l.decode()  
                q.put(newoutput)  

    def Running_Snakemake(self, event):
        try :             
            self.text1.SetLabelText("Started")
            self.progress_bar.Show()
            inputNf1 = 'snakemake -s ' + self.workspace + '/workflow/SnakeFile --cores ' + self.listbook.thread +' --use-conda'
            
            try :
                dir = '.snakemake/locks'
                dir2= '.snakemake/incomplete'
                filelist = glob.glob(os.path.join(dir, "*"))
                filelist2 = glob.glob(os.path.join(dir2, "*"))
                for f in filelist:
                    os.remove(f)
                for f in filelist2:
                    os.remove(f)

            except FileNotFoundError:
                pass

            self.listbook.fill_config_snakemake()

            self.process_sk = subprocess.Popen(inputNf1, shell=True,
                                            stdout = subprocess.PIPE,
                                            stderr = subprocess.STDOUT)
            
            pa_q = queue.Queue()
        
            pa_t = threading.Thread(target=self.read_output, args=(self.process_sk.stdout, pa_q))
            pa_t.daemon = True
            pa_t.start()

            while True: 
                self.progress_bar.Pulse()
                wx.GetApp().Yield()
                self.process_sk.poll()
                try:
                    l = pa_q.get(False)
                    self.listbook.TextCtrlRunSK.AppendText("\n" + str(l)) 
                except queue.Empty:
                    pass

                if self.process_sk.returncode is not None or self.process_sk.stdout== '' :
                    break
            
            self.text1.SetLabelText("Finished") 
            self.progress_bar.Hide()               

        except ValueError :
            dlg = ExceptionDialog("Please fill in the fields correctly !")
            dlg.ShowModal()
            if dlg.ShowModal == wx.ID_OK or wx.ID_CANCEL:
                dlg.Destroy()  

        event.Skip()



def main():
    app = wx.App(False)  # False pour ne pas initialiser un redirector d'affichage (console)
    ex = MyFrame(None, 'Esperanto')  # Créer une instance de la classe MyFrame
    ex.Show()  # Afficher la fenêtre principale
    app.MainLoop()  # Démarrer la boucle principale de l'application

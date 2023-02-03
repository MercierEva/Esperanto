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

path_to_home = os.path.expanduser('~')
path_to_local = os.path.join(path_to_home + '/.local')
path_to_gtk = os.path.join(path_to_local +'/share')

access_rights = 0o755

try:
    os.mkdir(path_to_local, access_rights)
    os.mkdir(path_to_gtk, access_rights)
except FileExistsError:
    pass

class Listbook(wx.Listbook):

    def __init__(self, parent):
        wx.Listbook.__init__(self, parent,  id=wx.ID_ANY, style=wx.EXPAND , size = (100,250))

        # panel pour le pipeline verticale
        self.panelpage1 = wx.Panel(self, wx.ID_ANY)
        self.panelpage1.SetBackgroundColour('#93c2cf')

        self.panelpage1bis = wx.Panel(self, wx.ID_ANY)
        self.panelpage1bis.SetBackgroundColour('#93c2cf')

        self.panelpage2 = wx.Panel(self, wx.ID_ANY)
        self.panelpage2.SetBackgroundColour('#93c2cf')

        il = wx.ImageList(30, 30)

        self.workspace = os.getcwd()

        img1 = wx.Image(self.workspace + "/images/clipart1544068.png")
        img1.Rescale(30, 30)
        btm1 = img1.ConvertToBitmap()
        il.Add(btm1)

        img1bis = wx.Image(self.workspace + "/images/clipart37384.png")
        img1bis.Rescale(30, 30)
        btm1bis = img1bis.ConvertToBitmap()
        il.Add(btm1bis)

        img2 = wx.Image(self.workspace + "/images/clipart2114274.png")
        img2.Rescale(30, 30)
        btm2 = img2.ConvertToBitmap()
        il.Add(btm2)

        self.AssignImageList(il)

        # Remplissage Panel page 1 !
        self.InsertPage(0, self.panelpage1, "Quality Control", imageId=0)
        
        vboxPanel1 = wx.BoxSizer(wx.VERTICAL)

        txttitle = ' Data length analysis '
        st2 = wx.StaticText(self.panelpage1, label=txttitle, style=wx.ALIGN_LEFT)
        self.fonttitle = wx.Font(13, family=wx.DECORATIVE, style=wx.NORMAL,  weight=wx.LIGHT)
        st2.SetFont(self.fonttitle)

        bar1 = wx.StaticLine(self.panelpage1, -1)

        vboxPanel1.Add(st2, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 15)
        vboxPanel1.Add(bar1, 0, wx.EXPAND | wx.BOTTOM | wx.RIGHT, 5)

        hboxP1_1 = wx.BoxSizer(wx.HORIZONTAL)  

        self.buttonload1 = wx.Button(self.panelpage1, label="Click to select archives of input fastq files :")
        self.buttonload1.SetBackgroundColour('#9EA3FF')

        hboxP1_1.Add(self.buttonload1, -1, wx.EXPAND |  wx.ALL, 5)
        vboxPanel1.Add(hboxP1_1, 0,  wx.EXPAND | wx.ALL, 5)

        self.listarr0 = MyListCtrl(self.panelpage1, -1)
        vboxPanel1.Add(self.listarr0, 1, wx.EXPAND | wx.ALL, border=5)
        dt = FileDrop(self.listarr0)
        self.listarr0.SetDropTarget(dt)    

        bar2 = wx.StaticLine(self.panelpage1, -1)
        vboxPanel1.Add(bar2, 0, wx.EXPAND | wx.RIGHT | wx.BOTTOM, 5)

        hboxP1_2 = wx.BoxSizer(wx.HORIZONTAL)
        txt_interval="max size of data : "
        stinterval = wx.StaticText(self.panelpage1, wx.ID_ANY, label=txt_interval)
        self.ctrl_param_interval = wx.TextCtrl(self.panelpage1, wx.ID_ANY,size=((50,17)))
        self.ctrl_param_interval.SetValue('6000')
        hboxP1_2.Add(stinterval, 0, wx.LEFT, 20)
        hboxP1_2.Add(self.ctrl_param_interval, 0, wx.LEFT, 5)

        self.interval= self.ctrl_param_interval.GetValue()

        self.buttonexe1 = wx.Button(self.panelpage1, label="Show histograms of lengths of data loading", size=(400, 25))
        self.buttonexe1.SetBackgroundColour('#f5d45e')
        
        hboxP1_2.Add(self.buttonexe1, 1, wx.LEFT, 20)

        vboxPanel1.Add(hboxP1_2)

        self.panelpage1.SetSizer(vboxPanel1)     
        self.panelpage1.Bind(wx.EVT_BUTTON, lambda event: self.IncludeFiles(event, self.listarr0), self.buttonload1)
        self.panelpage1.Bind(wx.EVT_TEXT, self.Change_name_P1)
        self.panelpage1.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, lambda event: self.OnRightDown(event, self.listarr0))


        #Remplissage Panel 1bis

        self.InsertPage(1, self.panelpage1bis, "Histograms Viewer", imageId=1)

        vboxpanel1bis = wx.BoxSizer(wx.VERTICAL)

        self.TextCtrlRunHisto = wx.TextCtrl(self.panelpage1bis, style= wx.TE_MULTILINE | wx.TE_RICH)
        self.TextCtrlRunHisto.SetInsertionPoint(0)
        self.TextCtrlRunHisto.SetBackgroundColour(wx.BLACK)
        self.TextCtrlRunHisto.SetForegroundColour(wx.WHITE)
        self.TextCtrlRunHisto.WriteText('Histograms of lengths: \n')
        vboxpanel1bis.Add(self.TextCtrlRunHisto, 1, wx.EXPAND | wx.ALL, 5)

        self.panelpage1bis.SetSizer(vboxpanel1bis)       

        # Remplissage Panel page 2
        self.InsertPage(2, self.panelpage2, "Treatment", imageId=2)
        vboxPanel2 = wx.BoxSizer(wx.VERTICAL)
        hboxssts = wx.BoxSizer(wx.HORIZONTAL)

        stt = wx.StaticBox(self.panelpage2, 0, "Step 1")
        sttSizer = wx.StaticBoxSizer(stt, wx.VERTICAL)
        sttBox = wx.BoxSizer(wx.VERTICAL)

        hboxstt1 = wx.BoxSizer(wx.HORIZONTAL)
        txt1 = 'Click to select archives of input fastq files :'
        st1 = wx.StaticText(self.panelpage2, label=txt1)
        loadButton = wx.Button(self.panelpage2, wx.ID_ANY, label='fastq.gz', size=((100,30)))
        loadButton.SetBackgroundColour('#0fb471')
        hboxstt1.Add(st1, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP,  border=15)
        hboxstt1.Add(loadButton, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP, 10 )
        sttBox.Add(hboxstt1, 0, wx.EXPAND | wx.ALL , border=10)
        
        hboxstt2 = wx.BoxSizer(wx.HORIZONTAL)
        stname_folder = wx.StaticText(self.panelpage2, wx.ID_ANY, label='Name of working folder : ')
        hboxstt2.Add(stname_folder, 0, wx.EXPAND | wx.ALIGN_LEFT | wx.LEFT|wx.TOP|wx.BOTTOM, border=15)   
        self.ctrl_name_folder = wx.TextCtrl(self.panelpage2, wx.ID_ANY, size=((180,30)))
        timestr = time.strftime("%Y%m%d") 
        self.ctrl_name_folder.SetValue('WorkingSpace_' + timestr)
        hboxstt2.Add(self.ctrl_name_folder, 0, wx.LEFT | wx.TOP , 10)
        self.killbutton=wx.Button(self.panelpage2, wx.ID_ANY, label="KILL Process", size=((110,30)))
        self.killbutton.SetBackgroundColour('#b40f52')
        hboxstt2.Add(self.killbutton, 0, wx.LEFT | wx.TOP, 10)

        sttBox.Add(hboxstt2, 0, wx.EXPAND | wx.ALL , border=10)
    
        sttSizer.Add(sttBox, 0, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL)

        stt2 = wx.StaticBox(self.panelpage2, 0, "Step 2")
        stt2Sizer = wx.StaticBoxSizer(stt2, wx.VERTICAL)
        stt2Box1 = wx.BoxSizer(wx.VERTICAL)

        stt2boxlength = wx.BoxSizer(wx.HORIZONTAL)      
        txt_minl = wx.StaticText(self.panelpage2, label="Minimum length : ", style = wx.ALIGN_CENTRE)
        self.ctrl_minl = wx.TextCtrl(self.panelpage2, wx.ID_ANY, size=((50,20)))
        txt_maxl = wx.StaticText(self.panelpage2, label="Maximum length : ",style = wx.ALIGN_CENTRE)
        self.ctrl_maxl = wx.TextCtrl(self.panelpage2, wx.ID_ANY, size=((50,20)))

        stt2boxlength.Add(txt_minl, 0, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_CENTER_HORIZONTAL, 5)
        stt2boxlength.Add(self.ctrl_minl, 0,wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL , 5)
        stt2boxlength.AddStretchSpacer()
        stt2boxlength.Add(txt_maxl,  0,wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL| wx.ALIGN_CENTER_HORIZONTAL, 5)
        stt2boxlength.Add(self.ctrl_maxl, 0,wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL, 5)
        stt2boxlength.AddStretchSpacer()


        stt2boxquality = wx.BoxSizer(wx.HORIZONTAL)
        txt_type = wx.StaticText(self.panelpage2, label="Read Type : ", style=wx.ALIGN_CENTRE)
        self.ctrl_type = wx.ComboBox(self.panelpage2, wx.ID_ANY, style=wx.CB_READONLY, choices=["1D", "2D", "1D2"], size=((60,20)))
        self.ctrl_type.SetValue("1D")

        txt_minimum_content= wx.StaticText(self.panelpage2, label="Min content (filtration) :")
        self.ctrl_min_cont = wx.TextCtrl(self.panelpage2, wx.ID_ANY, size=((40,20)))
        self.ctrl_min_cont.SetValue('75')
        txt_threading = wx.StaticText(self.panelpage2, label="Threads :",style = wx.ALIGN_CENTRE)
        self.ctrl_thread = wx.TextCtrl(self.panelpage2, wx.ID_ANY, size=((30,20)))
        
        stt2boxquality.Add(txt_type, 0,wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL| wx.ALIGN_CENTER_HORIZONTAL, 5)
        stt2boxquality.Add(self.ctrl_type, 0, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL, 5)
        stt2boxquality.Add(txt_minimum_content, 0,wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL| wx.ALIGN_CENTER_HORIZONTAL, 5)
        stt2boxquality.Add(self.ctrl_min_cont, 0,wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL, 5)
        stt2boxquality.Add(txt_threading, 0,wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL| wx.ALIGN_CENTER_HORIZONTAL, 5)
        stt2boxquality.Add(self.ctrl_thread, 0, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL, 5)


        stt2boxamorceR = wx.BoxSizer(wx.HORIZONTAL)
        txt_amorceR = wx.StaticText(self.panelpage2, label="To give the reverse primer (UPPER) : ", style = wx.ALIGN_CENTRE)
        self.ctrl_amorce = wx.TextCtrl(self.panelpage2, wx.ID_ANY, size=((250, 20)))
        self.ctrl_amorce.SetValue("")
        stt2boxamorceR.Add(txt_amorceR,  0, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL| wx.ALIGN_CENTER_HORIZONTAL, 5)
        stt2boxamorceR.Add(self.ctrl_amorce, 0, wx.ALL|wx.EXPAND, 5)

        stt2Box1.Add(stt2boxlength ,0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.EXPAND, border=5)
        stt2Box1.Add(stt2boxquality,0,wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.EXPAND, border=5)
        stt2Box1.Add(stt2boxamorceR, 0, wx.ALL | wx.EXPAND, border=5)

        stt2Sizer.Add(stt2Box1, 1, wx.EXPAND | wx.ALL, border = 5)

        hboxssts.Add(sttSizer, 1, wx.EXPAND | wx.ALL, 5)
        hboxssts.Add(stt2Sizer,1, wx.EXPAND | wx.ALL, 5)

        vboxPanel2.Add(hboxssts, 1, wx.EXPAND | wx.ALL, 5)
        self.buttonRun = wx.Button(self.panelpage2, wx.ID_ANY, label="Run Pipeline")
        self.buttonRun.SetBackgroundColour('#00ff7f')
        vboxPanel2.Add(self.buttonRun, 0, wx.EXPAND | wx.ALL , 5)

        self.minl = self.ctrl_minl.GetValue()
        self.maxl = self.ctrl_maxl.GetValue()
        self.type = self.ctrl_type.GetValue()
        self.thread = self.ctrl_thread.GetValue()
        self.amorce = self.ctrl_amorce.GetValue()
        self.folder = self.ctrl_name_folder.GetValue()
        self.content =self.ctrl_min_cont.GetValue()

        self.listarr1 = MyListCtrl(self.panelpage2, -1)
        vboxPanel2.Add(self.listarr1, 2, wx.EXPAND | wx.ALL, border=5)
        dt = FileDrop(self.listarr1)
        self.listarr1.SetDropTarget(dt)

        self.TextCtrlRunSK = wx.TextCtrl(self.panelpage2, style= wx.TE_MULTILINE | wx.TE_RICH)
        self.TextCtrlRunSK.SetInsertionPoint(0)
        self.TextCtrlRunSK.SetBackgroundColour(wx.BLACK)
        self.TextCtrlRunSK.SetForegroundColour(wx.WHITE)
        self.TextCtrlRunSK.AppendText("\n Here, the progression of workflow ")
        vboxPanel2.Add(self.TextCtrlRunSK, 2, wx.EXPAND | wx.ALL, 5)


        hboxbuttonfinaux=wx.BoxSizer(wx.HORIZONTAL)

        etiq_button_otu = " Compute statistic results of workflow "
        self.buttonLoadOtuTable = wx.Button(self.panelpage2, -1, label=etiq_button_otu)
        self.buttonLoadOtuTable.SetBackgroundColour('#f5d45e')
        hboxbuttonfinaux.Add(self.buttonLoadOtuTable, wx.ID_ANY, wx.EXPAND | wx.ALL, 5)

        vboxPanel2.Add(hboxbuttonfinaux, 0, wx.EXPAND | wx.ALL , 5)

        self.panelpage2.SetSizer(vboxPanel2)         

        self.panelpage2.Bind(wx.EVT_BUTTON, lambda event: self.IncludeFiles(event, self.listarr1), loadButton)
        self.panelpage2.Bind(wx.EVT_TEXT, self.Change_name_P2)
        self.panelpage2.Bind(wx.EVT_BUTTON, self.Build_OTU_table, self.buttonLoadOtuTable)
        self.panelpage2.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, lambda event: self.OnRightDown(event, self.listarr1))
    
    def Build_OTU_table(self, event):
        try:
            timestr = time.strftime("%Y%m%d-%H%M%S")
            list_of_files = glob.glob(self.workspace + '/workflow/'+ self.folder +'/05_stats/*final.tsv')

            if len(list_of_files) < 2 :
                try :
                    for one_file in list_of_files :
                        shutil.copy(one_file,  self.workspace + '/workflow/' + self.folder + '/05_stats/StatisticReport' + timestr + '.tsv')
                except :
                    pass
            else :
                li = []

                for filename in list_of_files:

                    df = pd.read_csv(filename, index_col=None, header=0, sep='\t' )

                    li.append(df)

                frame = pd.concat(li, axis=0, ignore_index=True)

                frame.to_csv( self.workspace + '/workflow/'+ self.folder +'/05_stats/StatisticReport' + timestr + '.tsv', sep='\t', header= True)



            subprocess.call(['soffice',  self.workspace + '/workflow/'+ self.folder +'/05_stats/StatisticReport' + timestr + '.tsv'])



        except ValueError:
            dlg = ExceptionDialog("This step requires that the workflow has been completed correctly !")
            dlg.ShowModal()
            if dlg.ShowModal == wx.ID_OK or wx.ID_CANCEL:
                dlg.Destroy()
    
    def Change_name_P2(self, event):

        self.minl = self.ctrl_minl.GetValue()
        self.maxl = self.ctrl_maxl.GetValue()
        self.type = self.ctrl_type.GetValue()
        self.thread = self.ctrl_thread.GetValue()
        self.amorce = self.ctrl_amorce.GetValue()
        self.folder = self.ctrl_name_folder.GetValue()
        self.content =self.ctrl_min_cont.GetValue()

        event.Skip()
    
    def IncludeFiles(self, event, array):

        try : 
            self.files = self.OnOpenFiles()

            for file_nano in self.files :
                basename = os.path.basename(file_nano)
                (name, ext) = os.path.splitext(basename)
                (name_cut, extcomp) = os.path.splitext(name)
                ex = extcomp + ext
                size = os.path.getsize(file_nano)
                sec = os.path.getmtime(file_nano)
                array.InsertItem(array.c, name_cut)
                array.SetItem(array.c, 1, ex)
                array.SetItem(array.c, 2, str(size) + ' B')
                array.SetItem(array.c, 3, time.strftime('%Y-%m-%d %H:%M', time.localtime(sec)))
                array.SetItem(array.c, 4, file_nano)

                if (array.c % 2) == 0:
                    array.SetItemBackgroundColour(array.c, '#e6f1f5')

                array.c += 1

            event.Skip()

        except:
            pass
    
    def Change_name_P1(self, event):

        self.interval = self.ctrl_param_interval.GetValue()
        event.Skip()  

    def OnOpenFiles(self):

        # otherwise ask the user what new file to open
        with wx.FileDialog(self, "Open fastq.gz file", wildcard="FASTQ files (*.fastq.gz)|*.fastq.gz",
                        defaultFile="*.fastq.gz",
                        style = wx.FD_OPEN | wx.FD_MULTIPLE) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_OK:
                # Proceed loading the file chosen by the user
                self.pathname = fileDialog.GetPaths()
                return self.pathname

            else :
                fileDialog.Destroy()
                return  # the user changed their mind

    
    def OnRightDown(self, event, array):
        menu = wx.Menu()
        suppr_one = menu.Append(wx.ID_ANY, 'remove one')
        suppr_all = menu.Append(wx.ID_ANY, 'remove all')
        menu.Bind(wx.EVT_MENU, lambda event: self.RemoveItem(event, array), suppr_one)
        menu.Bind(wx.EVT_MENU, lambda event: self.RemoveItems(event, array), suppr_all)
        self.PopupMenu(menu)
        menu.Destroy()
        event.Skip()
    
    def RemoveItem(self, event, array):
        array.DeleteItem(array.GetFocusedItem())
        array.c -=1
        event.Skip()
    
    def RemoveItems(self, event, array):
        array.DeleteAllItems()
        array.c = 0
        event.Skip()        

    def fill_config_snakemake(self):

        path_to_config = "workflow/config_wf.yaml"
        i=0
        dico_yaml = {}
        while i < self.listarr1.GetItemCount():
            name_file = self.listarr1.GetItemText(item=i, col=0)
            extension = self.listarr1.GetItemText(item=i, col=1)
            pathfile = self.listarr1.GetItemText(item=i, col=4)
            dico_yaml[str(name_file)] = str(pathfile)  
            i += 1

        dict_file = {'samples': dico_yaml , 
                    'folder' : 'workflow/' +self.folder+'/',
                    'params' : {'filtration' : {'min_length' : self.minl,
                    'max_length' : self.maxl,  
                    'readtype' : self.type, 
                    'min_content' : int(self.content)}, 
                    'threading' : int(self.thread)-1, 
                    'amorce_Reverse' : self.amorce}}
        
        with open(path_to_config, 'w') as file:
            yaml.dump(dict_file, file, default_flow_style=False, sort_keys=False)

class MyListCtrl(wx.ListCtrl):

    def __init__(self, parent, id):
        wx.ListCtrl.__init__(self, parent, id, style=wx.LC_REPORT )
        self.c=0

        self.InsertColumn(0, 'Name')
        self.InsertColumn(1, 'Ext')
        self.InsertColumn(2, 'Size', wx.LIST_FORMAT_RIGHT)
        self.InsertColumn(3, 'Modified')
        self.InsertColumn(4, 'path')

        self.SetColumnWidth(0, 130)
        self.SetColumnWidth(1, 70)
        self.SetColumnWidth(2, 100)
        self.SetColumnWidth(3, 130)
        self.SetColumnWidth(4, 70)
    
    def updateText(self, text):

        basename = os.path.basename(text)
        (name, ext) = os.path.splitext(basename)
        (name_cut, extcomp) = os.path.splitext(name)
        ex = extcomp + ext
        size = os.path.getsize(text)
        sec = os.path.getmtime(text)
        self.InsertItem(self.c, name_cut)
        self.SetItem(self.c, 1, ex)
        self.SetItem(self.c, 2, str(size) + ' B')
        self.SetItem(self.c, 3, time.strftime('%Y-%m-%d %H:%M', time.localtime(sec)))
        self.SetItem(self.c, 4, text)

        if (self.c % 2) == 0:
            self.SetItemBackgroundColour(self.c, '#e6f1f5')
        
        self.c += 1
    

class FileDrop(wx.FileDropTarget):

    def __init__(self, array):

        wx.FileDropTarget.__init__(self)
        self.array = array  

    def OnDropFiles(self, x, y, filenames):

        for each_file in filenames:
            self.array.updateText(each_file) 
        return True

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

def main():
    app = wx.App()
    ex = MyFrame(None, 'Esperanto')
    ex.Show()
    app.MainLoop()



if __name__ == '__main__':
    main()

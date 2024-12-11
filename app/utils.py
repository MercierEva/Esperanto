import os
import time
import gzip
import shutil
import wx
from .listctrl import MyListCtrl
from .filedrop import FileDrop
import yaml

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
        hboxstt1.Add(loadButton, 0, wx.TOP, 10 )
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
    
        sttSizer.Add(sttBox, 0, wx.EXPAND | wx.ALL)

        stt2 = wx.StaticBox(self.panelpage2, 0, "Step 2")
        stt2Sizer = wx.StaticBoxSizer(stt2, wx.VERTICAL)
        stt2Box1 = wx.BoxSizer(wx.VERTICAL)

        stt2boxlength = wx.BoxSizer(wx.HORIZONTAL)      
        txt_minl = wx.StaticText(self.panelpage2, label="Minimum length : ", style = wx.ALIGN_CENTRE)
        self.ctrl_minl = wx.TextCtrl(self.panelpage2, wx.ID_ANY, size=((50,20)))
        txt_maxl = wx.StaticText(self.panelpage2, label="Maximum length : ",style = wx.ALIGN_CENTRE)
        self.ctrl_maxl = wx.TextCtrl(self.panelpage2, wx.ID_ANY, size=((50,20)))

        stt2boxlength.Add(txt_minl, 0, wx.ALL|wx.EXPAND, 5)
        stt2boxlength.Add(self.ctrl_minl, 0,wx.ALL|wx.EXPAND, 5)
        stt2boxlength.AddStretchSpacer()
        stt2boxlength.Add(txt_maxl,  0,wx.ALL|wx.EXPAND, 5)
        stt2boxlength.Add(self.ctrl_maxl, 0,wx.ALL|wx.EXPAND, 5)
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
        
        stt2boxquality.Add(txt_type, 0,wx.ALL|wx.EXPAND, 5)
        stt2boxquality.Add(self.ctrl_type, 0, wx.ALL|wx.EXPAND, 5)
        stt2boxquality.Add(txt_minimum_content, 0,wx.ALL|wx.EXPAND, 5)
        stt2boxquality.Add(self.ctrl_min_cont, 0,wx.ALL|wx.EXPAND, 5)
        stt2boxquality.Add(txt_threading, 0,wx.ALL|wx.EXPAND, 5)
        stt2boxquality.Add(self.ctrl_thread, 0, wx.ALL|wx.EXPAND, 5)


        stt2boxamorceR = wx.BoxSizer(wx.HORIZONTAL)
        txt_amorceR = wx.StaticText(self.panelpage2, label="To give the reverse primer (UPPER) : ", style = wx.ALIGN_CENTRE)
        self.ctrl_amorce = wx.TextCtrl(self.panelpage2, wx.ID_ANY, size=((250, 20)))
        self.ctrl_amorce.SetValue("")
        stt2boxamorceR.Add(txt_amorceR,  0, wx.ALL|wx.EXPAND, 5)
        stt2boxamorceR.Add(self.ctrl_amorce, 0, wx.ALL|wx.EXPAND, 5)

        stt2Box1.Add(stt2boxlength ,0, wx.ALL | wx.EXPAND, border=5)
        stt2Box1.Add(stt2boxquality,0,wx.ALL | wx.EXPAND, border=5)
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

            list_of_files = glob.glob(self.workspace + '/workflow/'+ self.folder +'/07_stats/*final.tsv')

            if len(list_of_files) < 2 :
                try :
                    for one_file in list_of_files :
                        shutil.copy(one_file,  self.workspace + '/workflow/' + self.folder + '/StatisticReport.tsv')

                except :
                    pass
            else :
                li = []

                for filename in list_of_files:

                    df = pd.read_csv(filename, index_col=None, header=0, sep='\t' )

                    li.append(df)

                frame = pd.concat(li, axis=0, ignore_index=True)

                frame.to_csv( self.workspace + '/workflow/'+ self.folder +'/StatisticReport.tsv', sep='\t', header= True)



            subprocess.call(['soffice',  self.workspace + '/workflow/'+ self.folder +'/StatisticReport.tsv'])




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
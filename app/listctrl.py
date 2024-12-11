import os
import time
import wx


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
    

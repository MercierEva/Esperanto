import wx
import os


class FileDrop(wx.FileDropTarget):

    def __init__(self, array):

        wx.FileDropTarget.__init__(self)
        self.array = array  

    def OnDropFiles(self, x, y, filenames):

        for each_file in filenames:
            self.array.updateText(each_file) 
        return True
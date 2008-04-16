#!/usr/bin/python

# File: eden.py
from random import randint
from os import system
from Tkinter import *
from tkSimpleDialog import *
from tkMessageBox import *
from FileDialog import *

# When debug = 1 the script will print debug info!
debug = 1

tools = {
    'addmaps':      ['Add two binary maps', [ 'inp', 'bin', 'bin', 'bin'] ],
    'apodfc':	    ['Apodize the FCALC', [ 'inp', 'fcalc' ] ],
    'apodfo':	    ['Apodize the FOBS', [ 'inp', 'fobs' ] ],
    'back':	    ['Back transform an FCALC', [ 'inp', 'fcalc' ] ],
    'bin2map':	    ['Convert a .bin file to a .map file', [ 'inp' ] ],
    'cadhkl':	    ['Combine FCALC and FOBS files', [ 'inp' ] ],
    'count':	    ['Count electrons at atomic positions', [ 'inp' ] ],
    'distance':     ['Compute distance between maps', [ 'inp' ] ],
    'dphase':	    ['Compute phase difference between FCALCs', [ 'inp' ] ],
    'expandfc':     ['Expand the FCALC structure file to P1', [ 'inp' ] ],
    'expandfo':     ['Expand the FOBS structure file to P1', [ 'inp' ] ],
    'forth':	    ['Forward transform a binary map', [ 'inp' ] ],
    'map2bin':	    ['Convert .map file to .bin file', [ 'inp' ] ],
    'maketar':	    ['Generate a target', [ 'inp' ] ],
    'multmaps':     ['Multiply 2 binary maps', [ 'inp' ] ],
    'perturbhkl':   ['Add a perturbation to FCALC, real & imag.', [ 'inp' ] ],
    'ranphase':     ['Randomly perturb FCALC phases', [ 'inp' ] ],
    'regrid':	    ['Generate XPLOR .map file from Eden binary map', [ 'inp' ] ],
    'shapes':	    ['Find 3-D shapes in the density', [ 'inp' ] ],
    'solve':	    ['Optimize the density with respect to the data', [ 'inp', 'bin' ] ],
    'sym':	    ['Report symmetry-related pdb entries', [ 'inp' ] ],
    'tohu':	    ['Generate FCALC from pdb coordinates', [ 'inp' ] ],
    'variance':     ['Do statistics on the density', [ 'inp' ] ]
}
#    'bin2view':     'bin2view is no longer a part of Eden',
#    'view2bin':     'view2bin is no longer a part of Eden',
#    'drho':         '"drho" has been replaced by "distance"',
#    'setupncs':     '"setupncs" (for NCS option) has been disabled',
#    'tetra':        '"tetra" has been disabled.  Use "tohu" + "back"'

help_extra = {
   'keywords': 'Show a list of keywords',
   'news':     'Show latest news',
   'switches': 'Show command line switches'
   }


fff = {
    'new':       'New input file',
    'edit':      'Edit existing input file'
    }

class Callback:
    def __init__(self, function, *args, **kwargs):
        self.default_args = args
        self.default_kwargs = kwargs
        self.function = function
        
    def __call__(self, *more_args, **more_kwargs):
        # Create an updated argument tuple and kw-args dictionary.
        more_args = self.default_args + more_args
        more_kwargs.update(self.default_kwargs)

        # And call the function with all this stuff.
        self.function(*more_args, **more_kwargs)

class EdenInp:
    def __init__(self,home,root):
        """
        EDEN input file handling
        """
        self.home = home
        self.root = root
        try:
            edit = os.environ['EDITOR']
        except ValueError:
            nedit = [ 'gedit', 'nedit', 'kedit' ]
            ppp = os.path.split()
            for nn in nedit:
                for pp in ppp:
                    str = pp + '/' + nn
                    if os.path.exists(str):
                        edit = str
                        if debug == 1:
                            print "found "+str+ " \n"
        if edit is "vi" or edit is "vim":
            self.editor = "xterm -e " + edit
        else:
            self.editor = edit

    def edit_inp(self):
        tkf     = LoadFileDialog(self.root)
        self.fn = tkf.go("","*.inp","","")
        system(self.editor + " " + self.fn)
            
    def new_inp(self):
        # Command line options
        default = self.home + '/python/default.inp'
        system("cp " + default + " /tmp")
        system(self.editor + " /tmp/default.inp")
        tkf     = SaveFileDialog(self.root)
        self.fn = tkf.go("","*.inp","","")
        system("cp /tmp/default.inp " + self.fn) 

class MyAbout(Toplevel):
    def __init__(self,master,title,home):
        Toplevel.__init__(self,master)
        self.title(title)
        
        self.geometry(newGeometry="500x400")
        www = Canvas(self,background="White",relief=SUNKEN,
                     borderwidth=2,width=480,height=360)
        imagefile = home + '/python/ieden.gif'
        photo2    = PhotoImage(file=imagefile)
        item      = www.create_image(0,0,anchor=NW,image=photo2)
        www.pack(side=TOP)
        but       = Button(self,font=("Helvetica",16,"bold"),foreground="Navy",
                           text="(c) 2004, The EDEN team",command=self.destroy)
        but.pack(side=BOTTOM)
        self.mainloop()
    
class Eden:
    def __init__(self):
        """
        Interactive EDEN
        """

        self.title = " Interactive EDEN "
        self.homedir = os.environ['EDENHOME']

    def make_submenu(self,menu,header,mydict,cb):
        anamenu = Menu(menu)
        menu.add_cascade(label=header, menu=anamenu)
        for k,v in mydict.items():
            anamenu.add_command(label=k+' - '+v[0],command=Callback(cb,k))
        
    def make_gui(self):
        self.root = Tk(className=self.title)
        self.root.geometry(newGeometry="500x400")
        self.www = Canvas(self.root,background="White",relief=SUNKEN,
                          borderwidth=2,width=480,height=360)
        n = str(randint(1,2))
        
        imagefile = self.homedir + '/python/tyr' + n + '.gif'
        photo = PhotoImage(file=imagefile)
        item = self.www.create_image(0,0,anchor=NW,image=photo)
        self.www.pack(side=TOP)
        
        Button(self.root,font=("Helvetica",16,"bold"),foreground="Navy",
               text="http://www.edencrystallography.org",
               command=self.browse).pack(side=BOTTOM)
        #
        # Create the main menu
        #
        menu = Menu(self.root)
        self.root.config(menu=menu)
        #
        # Make the file menu
        #
        filemenu = Menu(menu)
        menu.add_cascade(label="File", menu=filemenu)
        for k,v in fff.items():
            filemenu.add_command(label=v,command=Callback(self.file_cb,k))
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.root.destroy)
        #
        # Make the tools menu
        #
        self.make_submenu(menu,"Tools",tools,self.run_cb)
#
#       Make the help menu
#
        helpmenu = Menu(menu)
        menu.add_cascade(label="Help", menu=helpmenu)
        helpmenu.add_command(label="About...", command=self.help_cb)
#        helpmenu.add_command(label="License", command=self.lic_cb)
        helpmenu.add_separator()
        for k,v in help_extra.items():
            helpmenu.add_command(label=k,command=Callback(self.th_cb,k))
        helpmenu.add_separator()
        for k,v in tools.items():
            helpmenu.add_command(label=k,command=Callback(self.th_cb,k))
        self.root.mainloop()

    def browse(self):
        system("mozilla http://www.edencrystallography.org &")
        
    def run_cb(self,k):
        # Here we need to insert code to generate the correct input
        # Preferably context sensitive
        if debug == 1:
            print "called the run_cb with flag ",k
        # Gather the command line options
        myprog = "xterm -e eden -v " + k
        for ff in tools[k][1]:
            tkf = LoadFileDialog(self.root)
            fn  = tkf.go("","*."+ff,"","")
            if fn == None:
                return
            myprog = myprog + " " + fn
#        myprog = myprog
#        tkf.done()
        if debug == 1:
            print "command line: ",myprog
        system(myprog)
        

    def file_cb(self,k):
        self.inp = EdenInp(self.homedir,self.root)
        if k is 'edit':
            self.inp.edit_inp()
        else:
            if k is 'new':
                self.inp.new_inp()
        
    def help_cb(self):
        info = MyAbout(self.www," About Interactive EDEN ",self.homedir)
        
    def th_cb(self,k):
        if debug == 1:
            print "called the th_cb with flag ",k
        info = "\n======= Calling information for " + k + " =======\n"
        cpfile = self.homedir + '/help/' + k + '.inf'
        f=open(cpfile, 'r')
        info = info + f.read()
        f.close()
        info = info + "   ======= Background information for " + k + " =======\n"
        cpfile = self.homedir + '/help/' + k + '.hlp'
        f=open(cpfile, 'r')
        info = info + f.read()
        f.close()
        info = info + "   ======= End of information for " + k + " =======\n"
        mmmm = Toplevel(self.root)
        scrollbar = Scrollbar(mmmm)
        scrollbar.pack(side=RIGHT,fill=Y)
        ttt = Text(mmmm,yscrollcommand=scrollbar.set,)
        ttt.insert(END,info)
        ttt.config(state=DISABLED)
        ttt.pack(side=LEFT)
        scrollbar.config(command=ttt.yview)

    def lic_cb(self):
        cpfile = self.homedir + '/python/COPYING'
        f=open(cpfile, 'r')
        license = f.read()
        f.close()
        mmmm = Toplevel(self.root)
        scrollbar = Scrollbar(mmmm)
        scrollbar.pack(side=RIGHT,fill=Y)
        ttt = Text(mmmm,yscrollcommand=scrollbar.set,)
        ttt.insert(END,license)
        ttt.config(state=DISABLED)
        ttt.pack(side=LEFT)
        scrollbar.config(command=ttt.yview)

# Now it's time to rock'n'roll...
godoit = Eden()
godoit.make_gui()


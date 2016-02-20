from Tkinter import *
import tkMessageBox

root = Tk()



def onFig7():
    fig7 = Toplevel(root, width=800, height=400)
    fig7.transient(root)
    fig7.grab_set()
    root.wait_window(fig7)



menu = Frame(root)
menu.pack(side=LEFT)

win = Frame(root,width=300, height=300)
win.pack(side=LEFT)

fig7btn = Button(menu, text="Figure 7", width=10, command=onFig7)
fig7btn.pack()

def callback():
    print "called the callback!"

def close():
    if tkMessageBox.askokcancel("Quit", "Do you really wish to quit?"):
        root.destroy()



root.protocol("WM_DELETE_WINDOW", close)
mainloop()

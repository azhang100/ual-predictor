


"""
def draw(bloodArray):
  for color in bloodArray:
    text = "{:5.2f}".format(color).zfill(5)
    if (len(text) > 5):
      text = text[:-1]
    print("|"+text, end="")
  print("|")

class Visualizer:
  
  def __init__(self, bloodArrayLen, timeArrayLen):
    self.bloodArrayLen = bloodArrayLen
    self.timeArrayLen = timeArrayLen
    # self.currTimeArrayIndex = 0
    
  def draw(self, bloodArray):
    assert(len(bloodArray) == self.bloodArrayLen)
    for bloodCellIndex in range(0,self.bloodArrayLen):
      color = bloodArray[bloodCellIndex]
      text = "|"+str(int(color)).zfill(2)
      print(text, end="")
    print("|")
"""

"""

from tkinter import *

class Visualizer(self){
  def __init__(self, bloodArrayLen, timeArrayLen){
    self.bloodArrayLen = bloodArrayLen
    self.timeArrayLen = timeArrayLen
    self.currTimeArrayIndex = 0
    self.top = Tk()
    
    C = Canvas(self.top, bg="blue", height=250, width=300)
    coord = 10, 50, 240, 210
    arc = C.create_arc(coord, start=0, extent=150, fill="red")
    line = C.create_line(10,10,200,200,fill='white')
    C.pack()
    self.top.mainloop()
    
  def draw(self, bloodArray){
    assert(len(bloodArray) == self.bloodArrayLen)
    for bloodCellIndex in range(0,self.bloodArrayLen):
      x = self.currTimeArrayIndex
      y = bloodCellIndex
      color = bloodArray[bloodCellIndex]
      your_canvas_widget.create_line(x, y, x + 1, y, fill=color)
}
"""

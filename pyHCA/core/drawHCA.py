#!/usr/bin/env python
""" The drawHCA module is an ensemble of functions to draw HCA plot of protein
sequences (with or without domain annotation
"""

import os, sys, argparse, string
import Bio.SeqIO
import Bio
from pyHCA.core.seq_util import transform_seq
from pyHCA.core.ioHCA import read_annotation

__author__ = "Tristan Bitard-Feildel"
__licence__= "MIT"
__version__ = 0.1
__email__ = "tristan [you know what] bitard@feildel.fr"
__institute__ = "IMPMC"

def readfseq(fseq):
    """
    brief read fasta sequence file and get header and sequence 
    param fseq is the address of fasta file
    return header is the header of fasta file
    return nbaa is the number of amino acid in the sequence
    return seq is the sequence as a string with 4 spaces at the begining and 4 at the end
    """
    
    f = open(fseq)
    seq= ""
    header = ""
    for i in f:
        if i[0] == ">":
            header = i.strip()
        
        else:
            seq+= i.strip()
    nbaa = len(seq)
    return header, nbaa, seq


def writeSVGheader(nbaa, fout):
    """
    brief write header of svg file
    param fout is an object file open in writing
    """
    
    fout.write("""<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="%d" height="210" version="1.1"
xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">

 <!-- ECMAScript with each click -->
  <script type="application/ecmascript"> <![CDATA[
    function aa_click(evt) {
      var aa = evt.target;
      var currentColor = aa.getAttribute("style");
      if (currentColor == "fill:blue; fill-opacity:0.0")
        {aa.setAttribute("style", "fill:blue; fill-opacity:0.3");
        }
      if (currentColor == "fill:blue; fill-opacity:0.3")
        {aa.setAttribute("style", "fill:red; fill-opacity:0.3");
        }
      if (currentColor == "fill:red; fill-opacity:0.3")
        {aa.setAttribute("style", "fill:gray; fill-opacity:0.3");
        }
      if (currentColor == "fill:gray; fill-opacity:0.3")
        {aa.setAttribute("style", "fill:blue; fill-opacity:0.0");
        }
        
    }
  ]]> </script>



"""%(((nbaa)/20.0)*90+150))
    
def getSVGheader(nbaa, height=210):
    """
    brief write header of svg file
    param fout is an object file open in writing
    """
    
    return ("""<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="{}" height="{}" version="1.1"
xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">

 <!-- ECMAScript with each click -->
  <script type="application/ecmascript"> <![CDATA[
    function aa_click(evt) {{
      var aa = evt.target;
      var currentColor = aa.getAttribute("style");
      if (currentColor == "fill:blue; fill-opacity:0.0")
        {{aa.setAttribute("style", "fill:blue; fill-opacity:0.3");
        }}
      if (currentColor == "fill:blue; fill-opacity:0.3")
        {{aa.setAttribute("style", "fill:red; fill-opacity:0.3");
        }}
      if (currentColor == "fill:red; fill-opacity:0.3")
        {{aa.setAttribute("style", "fill:gray; fill-opacity:0.3");
        }}
      if (currentColor == "fill:gray; fill-opacity:0.3")
        {{aa.setAttribute("style", "fill:blue; fill-opacity:0.0");
        }}
        
    }}
  ]]> </script>


""".format((((nbaa)/20.0)*90+150), height))

def getborder(x, y, F, yoffset=0):
    """
    """
    
    side  =    [[[ 0.83334, 7.555554,], [ -4.166666, 4.222221]],\
                [[-4.166666, 4.222221], [-2.5, -1.333335]],\
                [[-2.5, -1.333335], [4.166667, -3.555557]],\
                [[0.83334, 7.555554], [7.5, 5.333332]],\
                [[7.5, 5.333332],  [9.166667, -0.222224]],\
                [[9.166667, -0.222224], [4.166667, -3.555557]]] 
    
    lside = []
    polygon = "<polyline points=\""
    for begin, end in side:
        xb, yb, = begin
        xe, ye = end
        
        nxb = (x+xb)*F
        nyb = (y+yb*-1)*F
        nxe = (x+xe)*F
        nye = (y+ye*-1)*F
        
        polygon += "%f,%f %f,%f "%(nxb, nyb+yoffset, nxe, nye+yoffset)
        lside.append([[nxb,nyb], [nxe,nye]])
    polygon += """ " style="fill:black; fill-opacity:0.2" />"""
    return lside, polygon


def getCoordBorder(x, y, F=1, bis = 0, yoffset=0):
    """
    brief get the size limit of an amino acide in a HCAplot
    param x is the coord of the amino acid letter
    param y is the coord of the amino acid letter
    param F is the zoom index
    return lside is a list of 6 cases with coordinate of each x and y line around the amino acid in HCA plot
    """
    side  =    [[[ 0.83334, 7.555554,], [ -4.166666, 4.222221]],\
                [[-4.166666, 4.222221], [-2.5, -1.333335]],\
                [[-2.5, -1.333335], [4.166667, -3.555557]],\
                [[0.83334, 7.555554], [7.5, 5.333332]],\
                [[7.5, 5.333332],  [9.166667, -0.222224]],\
                [[9.166667, -0.222224], [4.166667, -3.555557]]] 
    
    
    lside = []
    polygon = "<polyline points=\""
    for begin, end in side:
        xb, yb, = begin
        xe, ye = end
        
        nxb = F * (x + xb)
        nyb = F * (y + yb) 
        nxe = F * (x + xe)
        nye = F * (y + ye) 
        lside.append([[nxb,nyb], [nxe,nye]])
        
        polygon += "%f,%f %f,%f "%(nxb+80, -nyb+80+bis+yoffset, nxe+80, -nye+80+bis+yoffset)
    
    polygon += """ " visibility="visible" style="fill:black; fill-opacity:0.0" onclick="aa_click(evt)"/>"""
    return lside, polygon
    

def getsideborder(seq, hydrophobe, nbaa):
    """
    brief get for each hydrophobe the cluster limit in a table
    param seq is the sequence with 4 spaces at the begining and at the end
    param hydrophobe is the string containing the hydrophobe amino acids
    param nbaa is the number of the initial protein
    return seqside is a table where each case contains 6 booleans representing the cluster limit
    """
    
    # adapted from the Luc Canard Algorithm
    
    # initialize of seqside to contains limit of each AA
    seqside = []
    for i in range(len(seq)):
        seqside.append([False, False, False, False, False, False])
    
    # for all position
    for n in range(0, len(seq)):
        
        if n >= 4 and n+1 < len(seq):       
            """ SIDE 0 """
            if seq[n] in hydrophobe:
                if seq[n-1] in hydrophobe:
                    seqside[n][0] = False
                else:
                    if seq[n-4] in hydrophobe:
                        if seq[n-3] == 'P' or seq[n-2] == 'P' or seq[n-1] == 'P':
                            seqside[n][0] = True
                        else:
                            seqside[n][0] = False
                    else:
                        seqside[n][0] = True
            else:
                if seq[n-4] in hydrophobe:
                    if seq[n-1] in hydrophobe:
                        if seq[n-3] == 'P' or seq[n-2] == 'P':
                            seqside[n][0] = False
                        else:
                            seqside[n][0] = True
                    else:
                        seqside[n][0] = False       
                        
            """ SIDE 1 """
            if ( seq[ n] in hydrophobe):      
                if ( seq[ n - 4] in hydrophobe):
                    if seq[ n - 3] == 'P' or seq[ n - 2] == 'P' or seq[ n - 1] == 'P':
                        seqside[ n][1] = True
                    else:
                        seqside[ n][1]= False
                else:
                    if ( seq[ n - 3] in hydrophobe ):
                        if seq[ n - 2]  == 'P' or  seq[ n - 1] == 'P':
                            seqside[ n][1] = True
                        else:    
                            seqside[ n][1] = False
                    else:
                        seqside[ n][1] = True
            else:
                if ( seq[ n - 4] in hydrophobe):
                    if ( seq[ n - 3] in hydrophobe):
                        seqside[ n][1] = True
                    else:
                        seqside[ n][1] = False
                else:
                    seqside[ n][1] = False
            
            """ SIDE 2"""
            if ( seq[ n] in hydrophobe):
                if ( seq[ n + 1] in hydrophobe):
                    seqside[ n][2] = False
                else:
                    if ( seq[ n - 3] in hydrophobe):
                        if ( seq[ n - 2] == 'P' or seq[ n - 1] == 'P'):    
                            seqside[ n][2] = True
                        else:
                            seqside[ n][2] = False
                    else:    
                        seqside[ n][2] = True
            else:
                if ( seq[ n - 3] in hydrophobe):
                    if ( seq[ n + 1] in hydrophobe):
                        if (    seq[ n - 2] == 'P' or seq[ n - 1] == 'P'  or seq[ n]    == 'P'    ):
                            seqside[ n][2] = False
                        else:
                            seqside[ n][2] = True
                    else:
                        seqside[ n][2] = False
                else:
                    seqside[ n][2] = False    

    
        if n < (len(seq)-4):
            """ SIDE 3"""
            if ( seq[ n] in hydrophobe ):
                if ( seq[ n - 1] in hydrophobe):
                    seqside[ n][3] = False
                else:
                    if (seq[ n + 3] in hydrophobe):
                        if ( seq[ n + 1]  == 'P' or seq[ n + 2] == 'P'):
                            seqside[ n][3] = True
                        else:
                            seqside[ n][3] = False
                    else:
                        seqside[ n][3] = True
            else:
                if ( seq[ n - 1] in hydrophobe):
                    if ( seq[ n + 3] in hydrophobe):
                        if (    seq[ n]  == 'P' or seq[ n + 1]  == 'P' or  seq[ n + 2] == 'P'):
                            seqside[ n][3] = False
                        else:
                            seqside[ n][3] = True
                    else:
                        seqside[ n][3] = False
                
                else:
                    seqside[ n][3] = False

            """ SIDE 4"""
            if ( seq[ n] in hydrophobe):
                if ( seq[ n +3] in hydrophobe):
                    if (     seq[ n + 1]  == 'P'   or seq[ n + 2] == 'P'):
                        seqside[n][4] = True
                    else:
                            seqside[n][4] = False
                else:
                    if ( seq[ n + 4] in hydrophobe):
                        if (    seq[ n + 1] == 'P' or seq[ n + 2] == 'P' or seq[ n + 3] == 'P'    ):
                            seqside[n][4] = True
                        else:
                            seqside[n][4] = False
                    else:
                        seqside[n][4] = True
            else:
                if n==1:
                    pass
                    #print n, seq[n], seq[n+3], seq[n+4]
                if ( seq[ n + 3] in hydrophobe):
                    if ( seq[ n + 4] in hydrophobe):
                        
                        seqside[n][4] = True
                        if n == 1 :
                            #print seqside[n]
                            pass
                    else:
                        seqside[n][4] = False
                else:
                    seqside[n][4] = False
        
            """ SIDE 5"""
            if ( seq[ n] in hydrophobe):
                if ( seq[ n + 1] in hydrophobe):
                    seqside[n][5] = False
                else:
                    if ( seq[ n + 4] in hydrophobe):
                        if ( seq[ n + 1] == 'P' or  seq[ n + 2] == 'P' or seq[ n + 3] == 'P' ):
                            seqside[n][5] = True
                        else:
                            seqside[n][5] = False
                    else:
                        seqside[n][5] = True
        
            else:
                if ( seq[ n + 4] in hydrophobe):
                    if ( seq[ n + 1] in hydrophobe):
                        if (    seq[ n + 2] == 'P'  or  seq[ n + 3] == 'P'):
                            seqside[n][5] = False
                        else:
                            seqside[n][5] = True
                    else:
                        seqside[n][5] = False
                else:
                    seqside[n][5] = False
    
    """ I do not understand this part
    around = [0,0,0,0,0,0,0,0, 0, 0]
    for n in range(nbaa+4, nbaa+2*4):
        
        around[ 0] = (seq[ n] in hydrophobe) ^ (seq[ n - 4] in hydrophobe)
        around[ 1] = (seq[ n]  in hydrophobe) ^ (seq[ n - 3]  in hydrophobe)
        around[ 3] = (seq[ n]  in hydrophobe) ^ (seq[ n - 1]  in hydrophobe)

        for i in range(2, 6):
            seqside[ n][i] = False
        
        if around[ 0]  and  around[ 3]:
            seqside[ n][ 0] = True
        else:
            seqside[ n][ 0] = False
        if around[ 0] and  around[ 3]:
            seqside[ n][ 1] = True
        else:
            seqside[ n][ 1] = False

    for n in range(0, 4):
        
        around[ 5] = (seq[ n] in hydrophobe) ^ (seq[ n + 4] in hydrophobe)
        around[ 7] = (seq[ n]  in hydrophobe) ^ (seq[ n + 3]  in hydrophobe)
        around[ 8] = (seq[ n]  in hydrophobe) ^ (seq[ n + 1]  in hydrophobe)

        for i in range(0, 4):
            seqside[ n][i] = False
        
        if around[ 7]  and  around[ 8]:
            seqside[ n][ 4] = True
        else:
            seqside[ n][ 4] = False
        if around[ 5] and  around[ 8]:
            seqside[ n][ 5] = True
        else:
            seqside[ n][ 5] = False
    """
    return seqside


    
def drawside(side, dybis = 0, second=80, yoffset=0): #seqside, lside, n, seq):
    """
    brief draw line for cluster in HCAplot
    param side is list of coord x and y
    param dybis is for mirror
    """
    
    begin, end = side
    x1, y1 = begin
    x2, y2 = end
    
    x1 += 80
    x2 += 80
    
    y1 = -y1 +(dybis+second)
    y2 = -y2 + (dybis+second)
    
    return """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:blue;stroke:black;stroke-width:1;" />"""%(x1,y1+yoffset, x2,y2+yoffset)
        

def deplace(i, seq, coord, side, hydrophobe, ):
    """
    brief deplace or not the line cluster to avoid hole in line cluster
    param i is the iterator on seq
    param seq is the sequence with 4 spaces at the begining and at the end
    param coord is a list with n case containing x,y coordinate of amino acid in hca plot
    param side is the current side of line cluster
    param hydrophobe is a list containing list of aa hydrophobe
    return answer is the distance to deplace line cluster
    """
    
    """ Algorithm by Luc Canard transcription from C program """
    
    before = i - 1
    after = i + 1

    answer = 0.
    dy = 0.
    dymoins3 = 0.
    dyplus3 = 0.
    dyplus4 = 0.
    dymoins4 = 0.
    HEIGTH = 40
    
    if ( seq[i] not in hydrophobe):
        #print len(coord), i, len(seq)
        
        # attention ici modifie avant 3 etait 4 et 2 etait 3
        dymoins4= coord[i][1] - coord[ i - 3][1]
        dyplus3 = coord[i][1] - coord[ i + 2][1] 
        dymoins3= coord[i][1] - coord[ i - 2][1]
        dyplus4 = coord[i][1] - coord[ i + 3][1]  
    
    #for j in range(6):
        if dymoins3 < 0.:
            dymoins3 = -dymoins3
        if dyplus3  < 0.:
            dyplus3 = -dyplus3
        if dymoins4 < 0.:
            dymoins4 = -dymoins4
        if dyplus4  < 0.:
            dyplus4 = -dyplus4
    
        if ( ( side == 0) and seq[ before] in hydrophobe and ( dymoins4 > 15.) ):
            dy = coord[i][1] - coord[before][1]
            if ( dy > 15.):
                answer = -HEIGTH
    
        if ( ( side == 3) and seq[ before] in hydrophobe and ( dyplus3 > 15.) ):
            dy = coord[i][1] - coord[ before][1]
            if ( dy > 15.):
                answer = -HEIGTH
        
        if ( ( side == 2) and seq[ after] in hydrophobe and ( dyplus4 > 15.)):
            dy = coord[i][1] - coord[ after][1]
            if ( dy < -15.):
                answer = HEIGTH
        
        if ( ( side == 5) and seq[ after] in hydrophobe and ( dymoins3 > 15.)):
            dy = coord[i][1] - coord[ after][1]
            if ( dy < -15.):
                answer = HEIGTH
    
    return answer

def getCoord(seq):
    """
    brief get coordinate of the plan helix for each position of the sequence 
    param seq is the sequence with 4 spaces at the end and at the begining
    return coord is a list with coordinate [x, y]
    """
    
    F = 1 # factor zoom
    FX = 12
    FY = 40
    AROUND = 4
    NTOUR = 3.6
    BULGARIANCONSTANT = 11
    dx = FX / NTOUR
    dy = FY / NTOUR
    
    coord = []
    #int_coord = []
    for n in range(len(seq)):
        
        x = dx * (n - AROUND)                           # +80
        y = -dy * ((n-AROUND+BULGARIANCONSTANT)%NTOUR)  #    +80
        ys = dy * ((n-AROUND+BULGARIANCONSTANT)%NTOUR)  #    +80 - 40
        
        coord.append([x, y])
        #ix = int(round(x))
        #iy = int(round(y))
        #int_coord.append((ix, iy))
        
    return coord #, int_coord


def drawProline(x, y, yoffset=0):
    """
    brief draw a star for proline
    param x coordinate
    param y coordinate
    return polygone is a string to make svg star
    """
    
    displace = [[350,75],\
                [379,161],\
                [469,161],\
                [397,215],\
                [423,301],\
                [423,301],
                [350,250],\
                [277,301],\
                [303,215],
                [231,161],\
                [321,161]]
    
    polygon = "<polyline points=\""
    for cur_x, cur_y in displace:
        
        x1 = (cur_x-321)/22.0 + x
        y1 = (cur_y-250)/22.0 + y
        
        polygon += "%f,%f "%(x1, y1+yoffset)
    
    polygon += """ " style=" fill:black; fill-opacity:1" />"""
    
    return polygon

def drawThreonine(x, y, yoffset=0):
    """
    brief draw a rectangle for threonine hca plot
    param x coord of x
    param y coord of y
    """
    
    rectangle = """<rect x="%f" y="%f" width="4" height="5" fill="none" stroke="black" stroke-width="0.2"/>"""%(x, y-5+yoffset)
    
    return rectangle


def drawSerine(x, y, yoffset=0):
    """
    brief draw a rectangle for serine hca plot
    param x coord of x
    param y coord of y
    """
    
    rectangle = """<rect x="%f" y="%f" width="4" height="5" fill="none" stroke="black" stroke-width="0.2"/>"""%(x, y-5+yoffset)
    prectangle = """<rect x="%f" y="%f" width="1" height="1" fill="black" stroke="black" stroke-width="0.2"/>"""%(x+1.5, y-3+yoffset)
    return rectangle+prectangle


def drawGlycine(x, y, yoffset=0):
    """
    brief draw a losange for glycine hca plot
    param x coord of x
    param y coord of y
    """
    displace = [[2,0],\
                [4,2],\
                [2,4],\
                [0,2]]
    
    polygon = "<polyline points=\""
    for cur_x, cur_y in displace:
        
        x1 = cur_x*1.2 + x+1
        y1 = cur_y*1.2 + y-5
        
        polygon += "%f,%f "%(x1, y1+yoffset)
    
    polygon += """ " style=" fill:black; fill-opacity:1" />"""
    return polygon


def linkCluster(seq, coords, n, dx, dy, yoffset=0):
    """
    brief trace the link to inside cluster
    param seq is a list containing residus and 4 spaces at the beg and at the end
    param coords is a list length of protein + 4 spaces in beg and 4 in end, each case contains [x, y] coordinates
    param n is the number of residu +4
    param dx to centre the draw
    param dy to centre the draw
    param fout svg opened file
    """
    
    dxmos = 4.0
    dymos = -3.5
    dxmos2 = 4.0
    dymos2 = 6.0
    
    svg = ""
    
    if coords[n+1][1] < coords[n][1]:
        x1 = coords[n][0]+dx+dxmos
        y1 = coords[n][1]+dy+dymos
        x2 = coords[n+1][0]+dx+dxmos
        y2 = coords[n+1][1]+dy+dymos
        svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)

        if coords[n+2][1] < coords[n+1][1]:
            x1 = coords[ n + 1][0] + dx + dxmos
            y1 = coords[ n + 1][1] + dy + dymos - 10
            x2 = ( coords[ n][0] + dx + dxmos + coords[ n + 1][0] + dx + dxmos) * 0.5
            y2 = ( coords[ n][1] + dy + dymos + coords[ n + 1][1] + dy + dymos - 40.) * 0.5
            svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
        else:
            x1 = coords[ n + 1][0] + dx + dxmos
            y1 = coords[ n + 1][1] + dy + dymos + 40. - 10
            x2 = ( coords[ n][0] + dx + dxmos + coords[ n + 1][0] + dx + dxmos) * 0.5
            y2 = ( coords[ n][1] + dy + dymos + coords[ n + 1][1] + dy + dymos + 40.) * 0.5
            svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
        
        x1 = coords[ n][0] + dx + dxmos
        y1 = coords[ n][1] + dy + dymos + 10
        x2 = ( coords[ n][0] + dx + dxmos + coords[ n + 1][0] + dx + dxmos) * 0.5
        y2 = ( coords[ n][1] + dy + dymos + coords[ n + 1][1] + dy + dymos + 40.) * 0.5
        svg +=  """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
    else:
        x1 = coords[ n][0] + dx + dxmos2
        y1 = coords[ n][1] + dy + dymos2
        x2 = coords[ n + 1][0] + dx + dxmos2
        y2 = coords[ n + 1][1] + dy + dymos2 - 20.
        svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
        
        x1 = coords[ n][0] + dx + dxmos
        y1 = coords[ n][1] + dy + dymos
        x2 = ( coords[ n][0] + dx + dxmos + coords[ n + 1][0] + dx + dxmos ) * 0.5
        y2 = ( coords[ n][1] + dy + dymos + coords[ n + 1][1] + dy + dymos - 40.) * 0.5
        svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
        
        x1 = ( coords[ n][0] + dx + dxmos + coords[ n + 1][0] + dx + dxmos )  * 0.5
        y1 = ( coords[ n][1] + dy + dymos + coords[ n + 1][1] + dy + dymos + 40.) * 0.5
        x2 = coords[ n + 1][0] + dx + dxmos
        y2 = coords[ n + 1][1] + dy + dymos
        svg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x1,-y1+yoffset, x2,-y2+yoffset)
        
    return svg
        

def dosvg(seq, coord, seqside, hydrophobe, nbaa, b, F=1, pathout=None, int_coords=None, yoffset=0):
    """
    brief draw a svg hca plot
    param seq is a list containing sequence and 4 spaces at the begining and at the end
    param coord is a list  [x, y] amino acid coordinate for each position
    param seqside is a boolean list containing 6 cases by case. 6 faces True or False to draw line cluster
    param nbaa is the number of residue
    param b if the begining if we use cutting option
    """
    
    """
    To make figure
    """
    """
    f = open("/home/guilhem/ARTICLE_DOMAIN/data/yeats/amasbreaker.res")
    for i in f:
        if i[0] == "H":
            hy =map(int, i.split()[1:])
        if i[0] == "B":
            br =map(int, i.split()[1:])
        if i[0] == "I":
            ic =map(int, i.split()[1:])
    f.close()
    """
    F = 1 # factor zoom
    FX = 12
    FY = 40
    AROUND = 4
    NTOUR = 3.6
    BULGARIANCONSTANT = 11
    dx = FX / NTOUR
    dy = FY / NTOUR

    coordbis = []
    for i, j in coord:
        coordbis.append([i, j-40])
    
    #im = Image.new("RGB", (,))
    
    if int_coords != None:
        int_coords2, int_coords_bis = [], []
        max_x, max_y = 0, 0
        max_x_bis, max_y_bis = 0, 0
        min_x, min_y = 1e10, 1e10
        min_x_bis, min_y_bis = 1e10, 1e10
        for x,y in int_coords:
            if max_x < x: 
                max_x = x
                max_x_bis = x
            if min_x > x: 
                min_x = x
                min_x_bis = x
            if max_y < y: 
                max_y = y
            if min_y > y: 
                min_y = y
            if max_y_bis < y - 40: 
                max_y_bis = y - 40
            if min_y_bis > y - 40: 
                min_y_bis = y - 40
            #print min_x, max_x, min_y, max_y
            #print x, y, y-40
            int_coords_bis.append((x/10, (y-40)/10))
            int_coords2.append((x/10, y/10))
        xlim = max(abs(max_x), abs(max_x_bis))
        int_coords2 = np.array(int_coords2)
        int_coords_bis= np.array(int_coords_bis)
        xlim = max(np.max(np.abs(int_coords2)), np.max(np.abs(int_coords_bis)))
        matrix = np.zeros( (9, xlim+1) )
        #print matrix.shape

    outsvg = ""
    
    for n in range(0, len(seq)):
        
        lside, polygon = getCoordBorder(coord[n][0], coord[n][1], yoffset=yoffset)
        lsides, polygon2 = getCoordBorder(coordbis[n][0], coordbis[n][1], yoffset=yoffset)
        
        x, y = coord[n]
        x += 80
        y = -y+80
        
        _, ys = coordbis[n]
        ys = -ys+80
        
        #print x, x2, y, ys
        
        # draw letter
        if seq[n] == "P":
            outsvg += drawProline(x, y, yoffset)
            outsvg += drawProline(x, ys, yoffset)
            if int_coords != None:
                matrix[ abs(int_coords2[n][1]), int_coords2[n][0] ] = -1
                matrix[ abs(int_coords_bis[n][1]), int_coords_bis[n][0] ] = -1
        elif seq[n] == "S":
            outsvg += drawSerine(x, y, yoffset)
            outsvg += drawSerine(x, ys, yoffset)
        elif seq[n] == "T":
            outsvg += drawThreonine(x, y, yoffset)
            outsvg += drawThreonine(x, ys, yoffset)
        elif seq[n] == "G":
            outsvg += drawGlycine(x, y, yoffset)
            outsvg += drawGlycine(x, ys, yoffset)
        elif seq[n] == "C":
            if int_coords != None:
                matrix[ abs(int_coords2[n][1]), int_coords2[n][0] ] = -2
                matrix[ abs(int_coords_bis[n][1]), int_coords_bis[n][0] ] = -2
            outsvg += """<text x="%f" y="%f" style="fill:black;font-size:8px;font-family:Times-Roman;stroke:black;stroke-width:0.4">%s</text>\n"""%(x*F, (y+yoffset)*F, seq[n])
            outsvg += """<text x="%f" y="%f" style="fill:black;font-size:8px;font-family:Times-Roman;stroke:black;stroke-width:0.4">%s</text>\n"""%(x*F, (ys+yoffset)*F, seq[n])
        else:
            outsvg += """<text x="%f" y="%f" style="fill:black;font-size:8px;font-family:Times-Roman">%s</text>\n"""%(x*F, (y+yoffset)*F, seq[n])
            outsvg += """<text x="%f" y="%f" style="fill:black;font-size:8px;font-family:Times-Roman">%s</text>\n"""%(x*F, (ys+yoffset)*F, seq[n])
            if seq[n] in hydrophobe and int_coords != None:
                matrix[ abs(int_coords2[n][1]), int_coords2[n][0] ] = 1
                matrix[ abs(int_coords_bis[n][1]), int_coords_bis[n][0] ] = 1
                
        # <ruler>
        if  n> 5 and ( n - 4 + 1) % ( 10) == 0 and n>=3:
            
            x =  coord[n][0] + 80 + 3.5#+80+dx
            dy1 = dy + coord[4][1] + 60 
            position = n - AROUND + 1 + b
            
            #print "AAAA",coord[n][0], dx
            #print "DIST", x-xp, x
            
            
            outsvg += """<text x="%f" y="%f" style="fill:black;font-size:8px;font-family:Times-Roman">%s</text>\n"""%(x*F, dy1-5+yoffset, position)
            outsvg += """<line x1="%f" y1="%f" x2="%f" y2="%f" style="fill:black;stroke:black;stroke-width:0.7;" />"""%(x, dy1+yoffset, x, dy1+5+yoffset)
        
        # cluster linker
        if seq[n] in hydrophobe and n+3 < len(seq):
            if (seq[n+1] not in hydrophobe and seq[n+2] in hydrophobe and
               seq[n+3] not in hydrophobe and seq[n+1] != "P" and
               seq[n-1] not in hydrophobe):
               outsvg += linkCluster(seq, coord, n, 80, -80, yoffset)
               outsvg += linkCluster(seq, coordbis, n, 80, -80, yoffset)
                    
        # each side of amino acid        
        for iteside, is_side in enumerate(seqside[n]):
            if is_side == False:
                continue
            
            if n+3 >= len(seq):
                continue
            dybis = deplace(n, seq, coord, iteside, hydrophobe)
            dybis2 = deplace(n, seq, coordbis, iteside, hydrophobe)
            
            outsvg += drawside(lside[iteside], -dybis, yoffset=yoffset)
            outsvg += drawside(lsides[iteside], -dybis2, yoffset=yoffset)
        
        """ for the figure
        if  not (n>3 and n<len(seq)-4):
            continue
        if n-4 in hy:
            polygon = polygon.replace("blue", "blue").replace("opacity:0.0", "opacity:0.2")
            polygon2 = polygon2.replace("blue", "blue").replace("opacity:0.0", "opacity:0.2")
            pass
        elif n-4 in ic:
            #polygon = polygon.replace("blue", "gray").replace("opacity:0.0", "opacity:0.2")
            #polygon2 = polygon2.replace("blue", "gray").replace("opacity:0.0", "opacity:0.2")
            pass
        else:
            polygon = polygon.replace("blue", "red").replace("opacity:0.0", "opacity:0.2")
            polygon2 = polygon2.replace("blue", "red").replace("opacity:0.0", "opacity:0.2")
            pass
        """
        outsvg += polygon
        outsvg += polygon2
    
    if int_coords != None:
        np.savez(pathout, mat=matrix)
    
    return outsvg

def domains2svg(start, stop, dom_name, status,  coord, yoffset=0):
    # ( n - 4 + 1) % ( 10) == 0 
    if start-1+4-1 > 0:
        x1 = coord[start+4-1][0] + 80 + 3.5#+80+dx
    else:
        x1 = coord[0][0] + 80 + 3.5#+80+dx
    if stop+4-1 < len(coord):
        x2 = coord[stop+4-1][0] + 80 + 3.5#+80+dx
    else :
        x2 = coord[-1][0] + 80 + 3.5
    size = x2 - x1
    
    name = """<text x="{}" y="{}" """.format(x1+10, yoffset+40)
    rect = """<rect x="{}" y="{}" """.format(x1, 25+yoffset)
    subrect = """<rect x="{}" y="{}" """.format(x1, 68+yoffset)
    if status == "!":
        name += """style="fill:black;font-size:12px;font-family:Times-Roman">"""
        rect += """style="fill:white; stroke:black; fill-opacity:1.0"  """
        subrect += """style="fill:white; stroke:black; fill-opacity:0.0"  """
    elif status == "o":
        name += """style="fill:black;font-size:12px;font-family:Times-Roman">"""
        rect += """style="fill:blue; stroke:black; fill-opacity:0.05"  """
        subrect += """style="fill:blue; stroke:black; fill-opacity:0.05"  """
    else:
        name += """style="fill:grey;font-size:12px;font-family:Times-Roman">"""
        rect += """style="fill:white; stroke:grey; stroke-dasharray=10,10; """
        subrect += """style="fill:white; stroke:grey; stroke-dasharray=10,10; """
        rect += """fill-opacity:1.0" """
        subrect += """fill-opacity:0.0" """
        
    name += """{}</text>\n""".format(dom_name)
    rect += """width="{}" height="20"/>\n\n""".format(size)
    subrect += """width="{}" height="100"/>\n\n""".format(size)
    return subrect+rect+name
    
    
def drawDomains(domains, coord, yoffset=0):
    """ From a list of domains (start, stop, names), draw a rectangle below 
    the HCA svg drawing
    
    Parameters
    ----------
    domains : list
        a list containing (start, stop, name) for each domains
    coord;
    yoffset:
    
    Return
    ------
    svg : string
    
    """
    return "\n".join([domains2svg(start, stop, name, status, coord, yoffset) 
                              for start, stop, name, nested, status in domains])

def draw_protnames(prot, yoffset=0):
    """ draw protein name in SVG
    """
    name = """<text x="0" y="{}" """.format(yoffset+25)
    name += """style="fill:black;font-size:12px;font-family:Times-Roman">"""
    name += """{}</text>\n""".format(prot)
    return name
    
def drawSelection(sequence, values, coord, yoffset=0):
    """ Draw pressure of selection in an svg format for each sequence
    
    Parameters
    ----------
    sequence
    values : list
        omega values
    coord;
    yoffset:
    
    Return
    ------
    """
    values_of_interest = []
    for i in range(len(sequence)):
        if sequence[i] != "-":
            values_of_interest.append(values[i])

    rect = ""
    for i, val in enumerate(values_of_interest):
        x = coord[i][0] + 80 + 3.5
        y = val
        rect += """<rect x="{}" y="{}" style="stroke:black; fill-opacity:0.0"  width="{}" height="20"/>\n\n""".format(x, y+yoffset, 1.0)

    return rect
    
def getSeqFromGi(gi):
    """
    brief get sequence from gi by fastacmd
    param gi is a gi number
    return temp_filename is the address of fasta file corresponding to the gi
    """
    FASTACMD = config.get("NCBI", "FASTACMD")
    NR = config.get("DATABASE", "NR")
    command = FASTACMD+" -d %s -s %s"%(NR, gi)
    #print command
    temp_fd, temp_filename = tempfile.mkstemp()
    process = subprocess.Popen(command, shell=True, stdout=temp_fd, 
                                                    stderr=subprocess.PIPE)
    process.wait()
    os.close(temp_fd)
    return temp_filename

def createHCAsvg(seq, nbaa, domains, b, yoffset=0):
    """
    brief draw the svg of hca
    param seq is the sequence with 4 spaces at both begining and end
    param nbaa is the true number of aa
    param output is the address of output file
    param b is the position to begin (for the ruler)
    """
    hydrophobe  = "YIMLFVW"
    seq = "    "+seq+"    "
    # n case with 6 sides in booleans
    seqside = getsideborder(seq, hydrophobe, nbaa)
    
    # get coord for amino acid in hca plot
    coord = getCoord(seq)
    #print len(coord), len(seq)
    #pathfig = output.split(".")[0]+".png"
    #fig, ax = plt.subplots()
    #for x,y in coord:
        #ax.set(aspect=1)
        #ax.plot(x,y-40,"o")
    #plt.savefig(pathfig)
    #plt.close()
    # draw the svg
    svg = dosvg(seq, coord, seqside, hydrophobe, nbaa, b, yoffset=yoffset)
    if domains:
        svg += drawDomains(domains, coord, yoffset)
    return svg
    
def drawHCA(prot, seq, nbaa, domains, output, b, yoffset=0):
    """
    brief draw the svg of hca
    param seq is the sequence with 4 spaces at both begining and end
    param nbaa is the true number of aa
    param output is the address of output file
    param b is the position to begin (for the ruler)
    """
    svgheader = getSVGheader(nbaa)
    namesvg = draw_protnames(prot, yoffset=yoffset)
    svg = createHCAsvg(seq, nbaa, domains, output, b, yoffset=yoffset)
    with open(output, "w") as outf:
        fout.write(svgheader)
        fout.write(namesvg)
        fout.write(svg)
        fout.write("</svg>")
    return int_coord


def make_svg(prot, prev_seq, annot, cnt):
    """ create svg for a given sequence
    """
    seq = transform_seq(prev_seq)
    nbaa = len(seq)
    
    b = 0
    svg = draw_protnames(prot, yoffset=cnt*230)
    svg += createHCAsvg(seq, nbaa, annot, b, yoffset=cnt*230)
    return svg, nbaa


def drawing(dfasta, annotation, pathout):
    """ draw multiple fasta sequences
    """
    svg = ""
    max_aa = 0
    cnt = 0
    for prot, prev_seq in dfasta.items():
        if type(prev_seq) == Bio.SeqRecord.SeqRecord:
            prev_seq = str(prev_seq.seq)
        cur_svg, nbaa = make_svg(prot, prev_seq, annotation.get(prot, []), cnt)
        svg += cur_svg
        if nbaa > max_aa:
            max_aa = nbaa
        cnt += 1
    
    # analys the new annotated domain, selective pressure from PAML
    #evolution_rate(pathnt, params.pathtree)
    svgheader = getSVGheader(max_aa, (cnt+1)*230)
    with open(pathout, "w") as fout:
        fout.write(svgheader)
        fout.write(svg)
        fout.write("</svg>")
    
def get_params():
    """ get command line ArgumentParser
    """
    parser = argparse.ArgumentParser(prog="{} {}".format(os.path.basename(sys.argv[0]), "draw"))
    parser.add_argument("-i", action="store", dest="fastafile", help="the fasta file", required=True)
    parser.add_argument("-d", action="store", dest="domain", help="[optionnal] provide domain annoation")
    parser.add_argument("-f", action="store", dest="domformat", help="the domain file format", choices=["pfam", "seghca"])
    parser.add_argument("-o", action="store", dest="svgfile", help="the svg file", required=True)
    params = parser.parse_args()
    return params

def main():
    # params 
    params = get_params()
    
    from pyHCA.core.ioHCA import read_multifasta

    # read fasta file
    dfasta = read_multifasta(params.fastafile)
    
    # compute hca annotation  
    annotation = {}
    if params.domain:
        annotation = read_annotation(params.domain, params.domformat)
        
    # draw
    drawing(dfasta, annotation, params.svgfile)
    
    sys.exit(0)
    
    
if __name__ == "__main__":
    main()



__version__ = "0.5"
__release__ = __version__  + '-dev' # comment out '-dev' before a release


from pyHCA.core.annotateHCA import main as annotate_main
from pyHCA.core.drawHCA import main as draw_main
from pyHCA.core.tremoloHCA  import main as tremolo_main


#def annotateHCA_main():
    #annotateHCA.main()
    
#if __name__ == "__main__" :
    #"""
    #Catching malfunctionning behavior
    #"""
    #try :
        #ext = annotateHCA_main( ) 
    #except SystemExit as e :
        #ext = e.code
    #except :
        #print("Unexpected error: {}".format(sys.exc_info()[0]))
        #traceback.print_exc()
        #ext = 1 
    #finally :    
        #sys.exit( ext )

__version__ = "0.5"
__release__ = __version__  + '-dev' # comment out '-dev' before a release


from .core.annotateHCA import main_segment as annotate_main
from .core.domainseq_segHCA import main_domainseq as domainseq_main
#from .core.annotateHCA import main_score as score_main
from .core.drawHCA import main as draw_main
from .core.tremoloHCA  import main as tremolo_main
from .core.domains_on_tree  import main as dom_ontree_main
from .core.HCA import HCA
from .core.disorderHCA import main as disorder_main

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

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__),"CommunFunctions"))
from ImageHandling import imread, imsave,SpatialImage
import numpy as np
from scipy import ndimage as nd
from cpp_wrapping import linearfilter, regionalext, connexe, watershed


def mars_segmentation(image_input, segmentation_output, sigma, h_min, sigma_ws, th=0, method=1, methods={1:'Classic',2:'Ace',3:'Hybridation'}, reconstructed_image=None,
    sigma_membrane=0.9, sensitivity=0.99, manual=False, manual_sigma=15, hard_thresholding=False, hard_threshold=1.0, sigma_TV=3.6, sigma_LF=0.9, sample=0.2):
    """ Perform the watershed segmentation of a given intensity image
    image_input : path to the input image
    segmentation_output : path to the final segmetented image
    sigma : value of the gaussian filter on the intensity image used for the h-minima operation (in voxels)
    h_min : value of the h-minima operation parameter (in bit)
    cleaning_needed : True if a removing of the cells touching the border of the image is needed
    sigma_ws : value of the gaussian filter on the intensity image used for the watershed (in voxel)
    th:Threshold for output cleaning

    method : int value which enables to switch the method for image_input preprocessing before mars segmentation process
            - 1 (default value): unchanged input image
            - 2 : 'Gace' reconstruction method applied on image_input. The output image is then sent to the mars segmentation process
            - 3 : 'Hybridation' method which processes a kind of merge between the original image image_input and its reconstruction with the 'Gace' reconstruction framework

    reconstructed_image : optional parameter for the path to be provided if the user wants to save the reconstructed image with the method.


    Ace method parameters:

        sigma_membrane=0.9 # parametre de rehaussement de membranes (en unites reelles, a priori 0.9 um est adapte a des images fusionnees de Patrick/Ralph/Aquila)

        # anisotropicHist /!\ etape critique
        sensitivity=0.99 # parametre de binarisation des membranes, etape critique /!\ en cas d'echec, privilegier une parametrisation "manuelle" de la fonction anisotropicHist via l'activation de l'option 'manual'

        manual=False     # par defaut, garder ce parametre a False. Si echec (ie seuils aberrants, se traduisant par une image binarisee des membranes de qualite tres mediocre), 
                         # le mettre a True et relancer les calculs sur l'image a tester. Si nouvel echec, jouer sur la valeur de manual_sigma... resultat non garanti
        manual_sigma=15  # parametre d'ajustement initial des histogrammes axiaux permettant le calcul de seuils de binarisation de l'image des membranes (parametre pris en compte seulement si manual = True). 
                         # Faire varier manual_sigma entre 5 et 25 en cas d'echec au premier essai. Resultat non garanti...

        hard_thresholding=False  # Si echec de la precedente methode de seuillage des membranes par calcul automatique de seuils directionnels, possibilite de choisir un seuillage dur global sur l'image en mettant cette option a True
        hard_threshold=1.0       # Si hard_thresholding est a True, seuillage des membranes rehaussees via ce seuil (valeur 1 : adaptee pour le time-point t001 d'Aquila par exemple)

        # TVmembrane
        sigma_TV=3.6     # parametre definissant l'echelle de vote (ie de propagation) des membranes par vote de tenseurs (en reelles), a priori a choisir entre 3 um (cellules etroites) et 4.5 um (gros gaps a combler)
        sigma_LF=0.9     # parametre de lissage de l'image reconstruite (en reelles), normalement la valeur par defaut = 0.9 um est satisfaisante
        sample=0.2       # parametre permettant d'optimiser la vitesse de calcul du vote de tenseurs (eviter de toucher)

    """
    program=mars_segmentation
    keep_reconstructed=True
    os.system('mkdir -p '+segmentation_output[:segmentation_output.rfind('/')])
    print 'Process Segmentation of '+image_input 
    if type(methods)==dict and methods.has_key(method):
        print "  Method '%s'"%methods[method]

    if image_input.split('.')[-1]=='tif' or image_input.split('.')[-1]=='tiff':
        print "%s: Warning: expectedly, input '%s' should not be in 'tif' format."%(program, image_input)
        print ' Convert in inr ' + image_input
        image_input_inr=''
        for n in image_input.split('.')[:-1]:
            image_input_inr+=n
        image_input_inr+='.inr'
        imsave(image_input_inr, imread(image_input))
        image_input=image_input_inr

    assert os.path.isfile(image_input), "%s: Error: file '%s' not found." %(program, image_input)
    
    assert reconstructed_image != image_input, "%s: the image_input and reconstructed_image parameters are both set to %s, which is forbidden. Exiting... "%(program,reconstructed_image)
    if method==1:
        if reconstructed_image==None:
            reconstructed_image=image_input
        else:
            imsave(reconstructed_image, imread(image_input))
            # Here keep_reconstructed = True on purpose, since we do not want to delete the original input data
    else:
        if reconstructed_image==None:
            reconstructed_image=segmentation_output.replace('.inr','_reconstructed_method_%s.inr'%(str(method)))
            keep_reconstructed=False
    
        if method==2:
            from ACE import GACE
            # Global Automated Cell Extractor
            out=GACE(image_input, path_output=reconstructed_image, sigma_membrane=sigma_membrane, sensitivity=sensitivity, manual=manual, manual_sigma=manual_sigma, hard_thresholding=hard_thresholding, hard_threshold=hard_threshold,
                     sigma_TV=sigma_TV, sigma_LF=sigma_LF, sample=sample, keep_all=False, keep_membrane=False, verbose=True)
        if method==3:
            print "%s: Error: the specified preprocessing method is not implemented yet. Please refer the project administrator for further information."%program
            quit()

    assert os.path.isfile(reconstructed_image), "%s: Error: file '%s' not found." %(program, reconstructed_image)

    #### Definition of paths to the different outputs ####
    path_gSigma=segmentation_output.replace('.inr','g' + str(sigma) + '.inr') # Path to the smoothed image for the local minima detection
    path_g_5=segmentation_output.replace('.inr','g'+str(sigma_ws)+'.inr') # Path to the smoothed image for the watershed
    path_rm=segmentation_output.replace('.inr','rm_s'+str(sigma)+'_h'+str(h_min)+'.inr') # Path to the regionalmax image
    path_cc=segmentation_output.replace('.inr','c_s'+str(sigma)+'_h'+str(h_min)+'.inr') # Path to the seeds image
    



    print 'Filter with sigma='+str(sigma)+' in ' + path_gSigma
    linearfilter(reconstructed_image, path_gSigma, sigma, realScale=True, type='deriche', verbose=False, lazy=True)
    
    print 'Filter with sigma='+str(sigma_ws)+' in ' + path_g_5
    linearfilter(reconstructed_image, path_g_5, sigma_ws, realScale=True, type='deriche', verbose=False, lazy=True)
    
    print 'Find Local minnima with h_min='+str(h_min)+' in ' + path_rm
    regionalext(path_gSigma, path_rm, h_min)
    
    print 'Find connex composant in ' + path_cc
    connexe(path_rm, path_cc, h_min)
    
    print 'Process watershed segmentation in ' + segmentation_output
    watershed(path_cc, path_g_5, segmentation_output)
    
    #Delete Temporary files
    cmd='rm ' + path_gSigma + ' ' + path_g_5 + ' ' + path_rm + ' ' + path_cc
    print cmd
    os.system(cmd)
    if not keep_reconstructed:
        cmd='rm ' + reconstructed_image
        print cmd
        os.system(cmd)
    print 'Segmentation done'


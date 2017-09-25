import numpy as np
import sys
import os
from scipy import ndimage as nd
import pickle as pkl
sys.path.append(os.path.join(os.path.dirname(__file__),"CommunFunctions"))
from ImageHandling import imread, imsave, SpatialImage
from cpp_wrapping import *





def random_number():
	import uuid
	return str(uuid.uuid4().fields[-1])[:5]+str(uuid.uuid4().fields[-1])[:5]



def light_LACE(parameters):
	"""
	LACE method with
	 - already transformed path_seg_0 (no need for path_fused_0 neither path_vector)
	 - already computed image of membranes (no need for path_fused_1 )
	Interface :
		(path_mask, label_of_interest, path_membrane_prefix, path_bin, rayon_dil, manual, manual_sigma, hard_thresholding, hard_threshold, sensitivity, verbose) = parameters
	 - path_mask : path to labelled image that defines the regions of interest (already transformed)
	 - label_of_interest : label of interest of path_mask image. The region of interest is then dilated with a ray of value rayon_dil.
	 - bbox : bounding box of the label of interest (before dilation) at the format [xmin, ymin, zmin, xmax, ymax, zmax].
	 - path_membrane_prefix : paths of maxima of enhanced membranes image (ext) and of associated angles (theta, phi) --> path_membrane_prefix+".[ext|theta|phi].inr" must be existing.
	 - path_bin : path to the binarised membranes result image
	 - rayon_dil : dilation ray for region of interest (see path_mask and label_of_interest)
	 - manual : parametre activant le mode manuel (ie parametrise) du seuil de binarisation des membranes si egal a True
	 - manual_sigma : parametre manuel d'initialisation pour le calcul du seuil de binarisation des membranes
	 - hard_thresholding : possibilite de choisir un seuillage dur global sur l'image en mettant cette option a True
	 - hard_threshold : si hard_thresholding est a True, seuillage des membranes rehaussees via ce seuil
	 - sensitivity : parametre de calcul des seuils anisotropiques selon un critere de sensibilite 
                                (true positive rate) : seuil = #(classe membrane>=seuil)/#(classe membrane)
	 - verbose : verbosite de la fonction
	light_LACE return value is path_bin.
	"""
	path_mask, label_of_interest, bbox, path_membrane_prefix, path_bin, rayon_dil, manual, manual_sigma, hard_thresholding, hard_threshold, sensitivity, verbose=parameters



	tmp_ID = random_number()	
	path_WORK = os.path.dirname(path_mask).rstrip(os.path.sep)+os.path.sep
	path_mask_dil = path_WORK+'mask_at_1_dil_'+tmp_ID+'.inr'

	path_ext_tmp = path_WORK + "tmp_"+tmp_ID+'.ext.inr'
	path_theta_tmp = path_WORK + "tmp_"+tmp_ID+'.theta.inr'
	path_phi_tmp = path_WORK + "tmp_"+tmp_ID+'.phi.inr'

	if not path_bin:
		path_bin = path_WORK + "tmp_"+tmp_ID+"."+str(label_of_interest)+'.inr'

	assert ( os.path.exists(path_membrane_prefix+".ext.inr")) and ( os.path.exists(path_membrane_prefix+".theta.inr")) and ( os.path.exists(path_membrane_prefix+".phi.inr"))

	if bbox:
		cropImage(path_mask, path_mask_dil, bbox, verbose=verbose )	
		seuillage(path_mask_dil, path_output=path_mask_dil,sb=label_of_interest, sh=label_of_interest, verbose=verbose )
	else:
		seuillage(path_mask, path_output=path_mask_dil,sb=label_of_interest, sh=label_of_interest, verbose=verbose )

	# Dilatation de la zone d'interet a l'instant t+1
	if True: # Calcul du rayon de dilatation en coordonnees reelles
		rayon_dil /= imread(path_mask).voxelsize[0]
		rayon_dil = int(rayon_dil+0.5)
	morpho(path_mask_dil, path_mask_dil, ' -dil -R '+str(rayon_dil), verbose=verbose)

	# Ici, path_mask_dil definit la zone d'interet dans laquelle LACE doit realiser la segmentation 

	if bbox:
		cropImage(path_membrane_prefix+".ext.inr", path_ext_tmp, bbox, verbose=verbose )

	# Binarisation des membranes
	if not hard_thresholding:
		# Seuillage anisotropique des membranes (parametre de sensitivite potentiellement critique)
		if bbox:
			cropImage(path_membrane_prefix+".theta.inr", path_theta_tmp, bbox, verbose=verbose )
			cropImage(path_membrane_prefix+".phi.inr", path_phi_tmp, bbox, verbose=verbose )
			anisotropicHist(path_input=path_ext_tmp, path_output=path_bin, path_mask=path_mask_dil, manual=manual, manual_sigma=manual_sigma, sensitivity=sensitivity, keepAll=False, verbose=verbose)
		else:
			anisotropicHist(path_input=path_membrane_prefix+".ext.inr", path_output=path_bin, path_mask=path_mask_dil, manual=manual, manual_sigma=manual_sigma, sensitivity=sensitivity, keepAll=False, verbose=verbose)
	else:
		if bbox:
			seuillage(path_input=path_ext_tmp, path_output=path_bin,sb=hard_threshold, verbose=verbose)
		else:
			seuillage(path_input=path_membrane_prefix+".ext.inr", path_output=path_bin,sb=hard_threshold, verbose=verbose)

	# Application du masque sur l'image binaire
	Logic(path_mask_dil, path_bin, path_bin, Mode='mask', verbose=verbose)


	if os.path.exists(path_mask_dil):
		cmd='rm ' + str(path_mask_dil)
		if verbose:
			print cmd
		os.system(cmd)
	if os.path.exists(path_ext_tmp):
		cmd='rm ' + str(path_ext_tmp)
		if verbose:
			print cmd
		os.system(cmd)
	if os.path.exists(path_theta_tmp):
		cmd='rm ' + str(path_theta_tmp)
		if verbose:
			print cmd
		os.system(cmd)
	if os.path.exists(path_phi_tmp):
		cmd='rm ' + str(path_phi_tmp)
		if verbose:
			print cmd
		os.system(cmd)
	
	return path_bin


def LACE(path_fused_0, path_fused_1, path_seg_0, label_of_interest, path_membrane_prefix=None, path_vector=None, path_bin=None, path_output=None, 
	rayon_dil=3.6, sigma_membrane=0.9, manual=False, manual_sigma=7, hard_thresholding=False, hard_threshold=1.0, sensitivity=0.99, sigma_TV=3.6, sigma_LF=0.9, sample=0.2, 
	keep_membrane=False, keep_hist=False, keep_vector=False, keep_all=False, short_LACE=False, verbose=False):
	'''
	LACE for Local Automated Cell Extractor
	

	# Paths d'entree 
	path_fused_0 : image fusionnee a l'instant t (celle dont on connait la segmentation) ; 
	path_fused_1 : image fusionnee a l'instant t+1 (celle qu'on souhaite reconstruire localement)
	path_seg_0 : image segmentee a l'instant t
	path_membrane_prefix+'.[ext|theta|phi].inr' (optionel) : paths des images de maxima (ext) de membranes rehaussees et de leurs angles associes (theta, phi), afin d'eviter de les recalculer
	path_vector (optionel) : path vers le champ de deformation calcule (via blockmatching) entre t (flo) et t+1 (ref) : T_flo<-ref
							 Si le path existe deja, le champ de deformation n'est pas recalcule, sinon il est calcule et conserve

	# Label d'interet de l'instant 0, definissant la zone d'etude propagee a l'instant 1
	label_of_interest : integer

	# Path de sortie 
	path_membrane_prefix+'.[ext|theta|phi].inr' (optionel) : paths de sauvegarde des images de maxima (ext) de membranes rehaussees et de leurs angles associes (theta, phi)
	path_vector (optionel) : cf paths d'entree
	path_bin (optionel) : path de sauvegarde de l'image des membranes binarisees (image envoyee en entree de l'etape de vote de tenseurs)
	path_output (optionel) : path de sauvegarde de l'image reconstruite de sortie (par defaut : None)

	# Mask parameters
	rayon_dil=3.6 (defaut, exprime en reelles) : rayon de dilatation de la zone d'etude propagee a partir de l'instant t vers l'instant t+1

	# Membrane reconstruction parameters
	sigma_membrane=0.9 (defaut, unites reelles, adapte a des images de telles que Patrick/Ralph/Aquila) : parametre de rehaussement des 
                                                                                   membranes avant binarisation de celles-ci
	sensitivity=0.99 (defaut) : parametre de calcul des seuils anisotropiques selon un critere de sensibilite 
                                (true positive rate) : seuil = #(classe membrane>=seuil)/#(classe membrane) 
	manual=False (defaut) : parametre activant le mode manuel (ie parametrise) du seuil de binarisation des membranes si egal a True
	manual_sigma=7 (defaut) : parametre manuel d'initialisation pour le calcul du seuil de binarisation des membranes

	hard_thresholding=False (defaut) : Si echec de la precedente methode de seuillage des membranes par calcul automatique de seuils directionnels,
                              possibilite de choisir un seuillage dur global sur l'image en mettant cette option a True
	hard_threshold=1.0 (defaut) : Si hard_thresholding est a True, seuillage des membranes rehaussees via ce seuil 
	                          (1.0 : adaptee pour le time-point t001 d'Aquila par exemple)

	sigma_TV=3.6 (defaut, reelles, adapte a des images de resolution 0.3um^3) : parametre de propagation des membranes 
                                            	                                par votes de tenseurs 
	sigma_LF=0.9 (defaut, reelles) : parametre de lissage gaussien de l'image des membranes reconstruite
	sample=0.2 (defaut) : echantillonne les votants de l'image initiale selon le coefficient (influe sur la vitesse de traitement)


	# version raccourcie de LACE (pour integration dans GLACE)
	short_LACE (default=False) : si True, ne procede pas au tensor voting et retourne l'image des membranes binarisees

	# Conserver ou non certaines images intermediaires
	keep_vector=False : conserver la transformation non lineaire T_t<-t+1 permettant de transformer l'image a t sur l'image a t+1 (se met a True automatiquement si path_vector est renseigne)
	keep_membrane=False : conserver toutes les images de l'etape de membranes (image d'extrema et d'angles) (se met a True automatiquement si path_membrane_prefix est renseigne)
	keep_hist=False : conserver le fichier d'histogrammes axiaux (nom de fichier se terminant par ".hist.txt")
	keep_all=False : conserver toutes les images intermediaires servant au calcul
	
	# Divers
	verbose=False : verbosite de la fonction, affiche par exemple les fichiers generes temporairement et ceux qui sont supprimes
	'''

	# Test existence des images d'entree
	assert(os.path.exists(path_fused_0) and os.path.exists(path_fused_1) and os.path.exists(path_seg_0))

	# Parametres de conservation des images generees
	if not keep_membrane:
		keep_membrane=keep_all
	if not keep_vector:
		keep_vector=keep_all
	if not keep_hist:
		keep_hist=keep_all

	# Definition des paths d'images intermediaires
	path_WORK=''
	if path_output:
		path_WORK = os.path.dirname(path_output).rstrip(os.path.sep)+os.path.sep
	else:
		if path_bin:
			path_WORK = os.path.dirname(path_bin).rstrip(os.path.sep)+os.path.sep
		else:
			path_WORK = os.path.dirname(path_fused_1).rstrip(os.path.sep)+os.path.sep
	if path_WORK==os.path.sep:
		path_WORK=os.path.curdir+os.path.sep

	tmp_ID = random_number()	
	if path_vector:
		keep_vector = True
	else:
		path_vector = path_WORK+'tmp_vectorfield_'+tmp_ID+'.inr'
	if path_membrane_prefix:
		keep_membrane = True
	else:
		path_membrane_prefix=path_WORK+'tmp_membrane_'+tmp_ID+''

	path_tmp = path_WORK+'tmp_threshold_'+tmp_ID+'.inr'
	path_mask = path_WORK+'mask_at_1_'+tmp_ID+'.inr.gz'
	path_mask_dil = path_WORK+'mask_at_1_dil_'+tmp_ID+'.inr.gz'
	path_affine_trsf = path_WORK+'tmp_affine_'+tmp_ID+'.txt'
	path_TV=path_WORK+'tmp_reconstructed_'+tmp_ID+'.inr'

	keep_tmp_bin=keep_all
	if not keep_tmp_bin:
		keep_tmp_bin=path_bin and (path_bin != path_membrane_prefix+'.bin.inr')

	# Verbose pour les fichiers
	if verbose:
		print "Temporary files:"
		print path_affine_trsf
		if not keep_vector:
			print path_vector
		print path_tmp
		print path_mask
		if not keep_membrane:
			if keep_tmp_bin:
				print path_membrane_prefix+".[ext|theta|phi].inr"
			else:
				print path_membrane_prefix+".[bin|ext|theta|phi].inr"
		else:
			if not keep_tmp_bin:
				print path_membrane_prefix+".bin.inr"
		if not keep_hist and not hard_thresholding:
			print path_membrane_prefix+".hist.txt"
		print path_TV
	if verbose and (keep_vector or path_output):
		print "Output files:"
		if keep_vector:
			print path_vector
		if keep_membrane:
			if keep_tmp_bin:
				print path_membrane_prefix+".[ext|theta|phi].inr"
			else:
				print path_membrane_prefix+".[bin|ext|theta|phi].inr"
		else:
			if keep_tmp_bin:
				print path_membrane_prefix+".bin.inr"
		if keep_hist and not hard_thresholding:
			print path_membrane_prefix+".hist.txt"
		if path_bin and path_bin != path_membrane_prefix+".bin.inr":
			print path_bin
		if path_output:
			print path_output


	### Path de sortie ###
	if not os.path.isdir(path_WORK):
		try:
			os.mkdir(path_WORK)
		except Exception :
			print "Unexpected error: unable to create working directory"

	### Calculs ###


	if not os.path.exists(path_vector):
		non_linear_registration(path_fused_0, path_fused_1, '/dev/null', path_affine_trsf, '/dev/null', path_vector, verbose=verbose)

	# Projection de la zone d'interet a l'instant t+1
	apply_trsf(path_seg_0, path_trsf=path_vector, path_output=path_mask, template=path_fused_1, nearest=True, verbose=verbose)


	if (not os.path.exists(path_membrane_prefix+".ext.inr")) or (not os.path.exists(path_membrane_prefix+".theta.inr")) or (not os.path.exists(path_membrane_prefix+".phi.inr")):
		# Extraction de la zone d'interet a l'instant t
		seuillage(path_mask, path_output=path_mask_dil,sb=label_of_interest, sh=label_of_interest, verbose=verbose )
		# Calcul du rayon de dilatation en coordonnees reelles
		rayon_dil_vox = rayon_dil / imread(path_mask).voxelsize[0]
		rayon_dil_vox = int(rayon_dil_vox+0.5)
		morpho(path_mask_dil, path_mask_dil, ' -dil -R '+str(rayon_dil_vox), verbose=verbose)
		# Renforcement local des membranes de l'image fusionnee a l'instant t+1 + extraction des maxima directionnels
		membrane_renforcement(path_fused_1, prefix_output=path_membrane_prefix, path_mask=path_mask_dil,  init=sigma_membrane, verbose=verbose)
		if path_mask_dil and os.path.exists(path_mask_dil):
			cmd='rm ' + path_mask_dil 
			if verbose:
				print cmd
			os.system(cmd)


	if not path_bin:
		path_bin=path_membrane_prefix+'.bin.inr'

	parameters=(path_mask, label_of_interest, None, path_membrane_prefix, path_bin, rayon_dil, manual, manual_sigma, hard_thresholding, hard_threshold, sensitivity, verbose)

	if verbose:
		print 'Running light_LACE(' + str(parameters) + ') ...'

	path_bin=light_LACE(parameters)

	reconstructed_image_1=None
	if short_LACE:
		if path_bin:
			reconstructed_image_1=path_bin
		else:
			reconstructed_image_1=path_membrane_prefix+'.bin.inr'
	else:
		# Tensor voting sur l'image des membranes binarisees
		TVmembrane(path_input=path_bin, path_output=path_TV, sample=sample, scale=sigma_TV, sigma_LF=sigma_LF, realScale=True, keepAll=False, verbose=verbose)

		# Copie vers l'image de sortie 
		if path_output and path_TV != path_output:
			copy(path_TV, path_output, verbose=verbose)

		reconstructed_image_1=path_TV

	# Suppression des images intermediaires
	files_to_rm = ""

	if not keep_all:
		if path_tmp and os.path.exists(path_tmp):
			files_to_rm += path_tmp + " "
		if path_mask and os.path.exists(path_mask):
			files_to_rm += path_mask + " "
		if path_affine_trsf and os.path.exists(path_affine_trsf):
			files_to_rm += path_affine_trsf + " "
		if path_vector and os.path.exists(path_vector) and not keep_vector:
			files_to_rm += path_vector + " "
		if path_TV and path_TV != path_output and os.path.exists(path_TV):
			files_to_rm += path_TV + " "

	if not keep_membrane:
		if path_membrane_prefix:
			if os.path.exists(path_membrane_prefix+'.ext.inr'):
				files_to_rm += path_membrane_prefix+'.ext.inr' + " "
			if os.path.exists(path_membrane_prefix+'.theta.inr'):
				files_to_rm += path_membrane_prefix+'.theta.inr' + " "
			if os.path.exists(path_membrane_prefix+'.phi.inr'):
				files_to_rm += path_membrane_prefix+'.phi.inr' + " "

	if not keep_tmp_bin and path_membrane_prefix and os.path.exists(path_membrane_prefix+'.bin.inr'):
		files_to_rm +=  path_membrane_prefix+'.bin.inr' + " "
		
	if not keep_hist and path_membrane_prefix and os.path.exists(path_membrane_prefix+'.hist.txt'):
		files_to_rm += path_membrane_prefix+'.hist.txt' + " "

	if files_to_rm:
		if verbose:
			print "Deleting temporary files: \n" + files_to_rm
		os.system("rm " + files_to_rm)

	# Fin de LACE
	return reconstructed_image_1



def GACE(path_input, binary_input=False, path_membrane_prefix=None, path_bin=None, path_output=None, 
		 sigma_membrane=0.9, manual=False, manual_sigma=7, sensitivity=0.99, hard_thresholding=False, hard_threshold=1.0, 
		 sigma_TV=3.6, sigma_LF=0.9, sample=0.2,
		 keep_membrane=False, keep_hist=False, keep_all=False, verbose=False):
	'''
	GACE for Global Automated Cell Extractor
	

	# Paths d'entree 
	path_input : image fusionnee OU image binaire que l'on souhaite reconstruire 
				 /!\ Si path_input est une image binaire utilisee avec l'option binary_input=True, veiller a ce que ces conditions soient verifiees :
						  - path_input = <prefix>.inr[.gz] ou <prefix>.<particule>.inr[.gz]
						  - il existe <prefix>.theta.inr et <prefix>.phi.inr qui correspondent aux images d'angles associees a l'orientation des membranes dans l'image binaire

	path_membrane_prefix+'.[ext|theta|phi].inr' (optionel) : paths des images de maxima (ext) de membranes rehaussees et de leurs angles associes (theta, phi), afin d'eviter de les recalculer (seulement si binary_input=False)

	# Path de sortie 
	path_membrane_prefix+'.[ext|theta|phi].inr' (optionel) : paths de sauvegarde des images de maxima (ext) de membranes rehaussees et de leurs angles associes (theta, phi) (seulement si binary_input=False)
	path_bin (optionel) : path de sauvegarde de l'image des membranes binarisees (image envoyee en entree de l'etape de vote de tenseurs)
	path_output (optionel) : path de sauvegarde de l'image reconstruite de sortie (par defaut : None)

	# Membrane reconstruction parameters
	sigma_membrane=0.9 (defaut, unites reelles, adapte a des images de telles que Patrick/Ralph/Aquila) : parametre de rehaussement des 
                                                                                   membranes avant binarisation de celles-ci

	sensitivity=0.99 (defaut) : parametre de calcul des seuils anisotropiques selon un critere de sensibilite 
                                (true positive rate) : seuil = #(classe membrane>=seuil)/#(classe membrane) 
	manual=False (defaut) : parametre activant le mode manuel (ie parametrise) du seuil de binarisation des membranes si egal a True
	manual_sigma=7 (defaut) : parametre manuel d'initialisation pour le calcul du seuil de binarisation des membranes
	
	hard_thresholding=False (defaut) : Si echec de la precedente methode de seuillage des membranes par calcul automatique de seuils directionnels,
                              possibilite de choisir un seuillage dur global sur l'image en mettant cette option a True
	hard_threshold=1.0 (defaut) : Si hard_thresholding est a True, seuillage des membranes rehaussees via ce seuil 
	                          (1.0 : adaptee pour le time-point t001 d'Aquila par exemple)

	sigma_TV=3.6 (defaut, reelles, adapte a des images de resolution 0.3um^3) : parametre de propagation des membranes 
                                                                              	par votes de tenseurs 
	sigma_LF=0.9 (defaut, reelles) : parametre de lissage gaussien de l'image des membranes reconstruite
	sample=0.2 (defaut) : echantillonne les votants de l'image initiale selon le coefficient (influe sur la vitesse de traitement)

	# Conserver ou non la transformation non lineaire
	keep_all=False : conserver toutes les images intermediaires servant au calcul
	keep_membrane=False : conserver toutes les images de l'etape de membranes (image d'extrema, d'extrema binarises et d'angles)
	
	# Divers
	verbose=False : verbosite de la fonction, affiche par exemple les fichiers generes temporairement et ceux qui sont supprimes
	'''

	# Test existence de l'image d'entree
	assert(os.path.exists(path_input))

	# Parametres de conservation des images generees
	if not keep_membrane:
		keep_membrane=keep_all

	# Definition des paths d'images intermediaires
	path_WORK=''
	if path_output:
		path_WORK = os.path.dirname(path_output).rstrip(os.path.sep)+os.path.sep
	else:
		path_WORK = os.path.dirname(path_input).rstrip(os.path.sep)+os.path.sep
	if path_WORK==os.path.sep:
		path_WORK=os.path.curdir+os.path.sep
	tmp_ID = random_number()	
	if path_membrane_prefix:
		keep_membrane = True
	else:
		path_membrane_prefix=path_WORK+'tmp_membrane_'+tmp_ID+''
	path_TV=path_WORK+'tmp_reconstructed_'+tmp_ID+'.inr'

	keep_tmp_bin=keep_all
	if not keep_tmp_bin:
		keep_tmp_bin=path_bin and (path_bin != path_membrane_prefix+'.bin.inr')


	if verbose:
		print "Temporary files:"
		if binary_input:
			print path_TV
		else:
			if keep_membrane:
				print path_membrane_prefix+".bin.inr"
			else:
				print path_membrane_prefix+".[bin|ext|theta|phi].inr"
			if not hard_thresholding:
				print path_membrane_prefix+".hist.txt"
			print path_TV
	if verbose and path_output:
		print "Output file:"
		if path_output:
			print path_output

	### Path de sortie ###
	if not os.path.isdir(path_WORK):
		try:
			os.mkdir(path_WORK)
		except Exception :
			print "Unexpected error: unable to create working directory"

	### Calculs ###
	
	path_TV_input=''
	if binary_input:
		# Image d'entree = image binaire
		path_TV_input=path_input
		keep_membrane=True
	else:
		if (not os.path.exists(path_membrane_prefix+".ext.inr")) or (not os.path.exists(path_membrane_prefix+".theta.inr")) or (not os.path.exists(path_membrane_prefix+".phi.inr")):
			# Renforcement des membranes de l'image fusionnee a l'instant t+1 + extraction des maxima directionnels
			membrane_renforcement(path_input, prefix_output=path_membrane_prefix, path_mask=None,  init=sigma_membrane, verbose=verbose)

		# Binarisation des membranes
		if not hard_thresholding:
			# Seuillage anisotropique des membranes (parametre de sensitivite potentiellement critique)
			anisotropicHist(path_input=path_membrane_prefix+".ext.inr", path_output=path_membrane_prefix+'.bin.inr', path_mask=None, manual=manual, manual_sigma=manual_sigma, sensitivity=sensitivity, keepAll=False, verbose=verbose)
		else:
			# Seuillage dur
			seuillage(path_input=path_membrane_prefix+".ext.inr", path_output=path_membrane_prefix+'.bin.inr',sb=hard_threshold, verbose=verbose)
		path_TV_input=path_membrane_prefix+".bin.inr"
		if path_bin and not os.path.exists(path_bin):
			# Copie de l'image binaire temporaire vers le path renseigne en parametre
			assert os.path.exists(path_membrane_prefix+'.bin.inr')
			copy(path_membrane_prefix+".bin.inr", path_bin, verbose=verbose)

	# Tensor voting sur l'image des membranes binarisees
	if verbose:
		print 'Processing Tensor Voting on image ' + path_TV_input + ' ...'
	assert(os.path.exists(path_TV_input))
	assert (path_TV_input.endswith('.inr.gz') or path_TV_input.endswith('.inr') )
	binary_file_decomp = path_TV_input.split('.')
	if binary_file_decomp[-1]=='gz':
		binary_file_decomp=binary_file_decomp[:-1]
	assert binary_file_decomp[-1]=='inr'
	binary_file_decomp=binary_file_decomp[:-1]
	if (not os.path.exists('.'.join(binary_file_decomp)+'.theta.inr')) or (not os.path.exists('.'.join(binary_file_decomp)+'.phi.inr')):
		assert (len(binary_file_decomp)>1 and os.path.exists('.'.join(binary_file_decomp[:-1])+'.theta.inr') and os.path.exists('.'.join(binary_file_decomp[:-1])+'.phi.inr')), "Error : unexpectedly, <prefix>.theta.inr and/or <prefix>.phi.inr not found for file " + path_TV_input + " before tensor voting step"

	TVmembrane(path_input=path_TV_input, path_output=path_TV, path_mask=None, scale=sigma_TV, sample=sample, sigma_LF=sigma_LF, realScale=True, keepAll=False, verbose=verbose)


	# Lecture de l'image reconstruite (image retournee par la fonction)
	reconstructed_image=imread(path_TV)

	# Copie vers l'image de sortie 
	if path_output and path_TV != path_output:
		copy(path_TV, path_output, verbose=verbose)

	# Suppression des images intermediaires
	files_to_rm = ""

	if not keep_all:
		if path_TV and path_TV != path_output and os.path.exists(path_TV):
			files_to_rm += path_TV + " "

	if not binary_input and not keep_membrane:
		if path_membrane_prefix:
			path_membrane_dir=os.path.dirname(path_membrane_prefix)
			files_to_rm += ''.join([path_membrane_dir+os.path.sep+i+" " for i in os.listdir(path_membrane_dir) if os.path.isfile(path_membrane_dir+os.path.sep+i) and (path_membrane_dir+os.path.sep+i).count(path_membrane_prefix)])


	if not keep_membrane:
		if path_membrane_prefix:
			if os.path.exists(path_membrane_prefix+'.ext.inr'):
				files_to_rm += path_membrane_prefix+'.ext.inr' + " "
			if os.path.exists(path_membrane_prefix+'.theta.inr'):
				files_to_rm += path_membrane_prefix+'.theta.inr' + " "
			if os.path.exists(path_membrane_prefix+'.phi.inr'):
				files_to_rm += path_membrane_prefix+'.phi.inr' + " "

	if not keep_tmp_bin and path_membrane_prefix and os.path.exists(path_membrane_prefix+'.bin.inr'):
		files_to_rm +=  path_membrane_prefix+'.bin.inr' + " "
		
	if not keep_hist and path_membrane_prefix and os.path.exists(path_membrane_prefix+'.hist.txt'):
		files_to_rm += path_membrane_prefix+'.hist.txt' + " "


	if files_to_rm:
		if verbose:
			print "Deleting temporary files: \n" + files_to_rm
		os.system("rm " + files_to_rm)

	# Fin de LACE
	return reconstructed_image


def GLACE(path_fused_0, path_fused_1, path_seg_0, labels_of_interest='all', background=[0,1], path_membrane_prefix=None, path_vector=None, path_bin=None, path_output=None, rayon_dil=3.6, 
	sigma_membrane=0.9, manual=False, manual_sigma=7, hard_thresholding=False, hard_threshold=1.0, sensitivity=0.99, sigma_TV=3.6, sigma_LF=0.9, sample=0.2, 
	keep_membrane=False, keep_vector=False, keep_all=False,  nb_proc=7, verbose=False):
	'''
	GLACE --- Grouped Local Automated Cell Extractor

	Summary of the method steps :
		# Step 1 : non-linear registration between path_fused_0 and path_fused_1

		# Step 2 : membrane reinforcement on global image

		# Step 3 : init and build the binary image (loop of short LACEs)

			# Step 3.1 : LACE on the label

			# Step 3.2 : adding to the main binary image

		# Step 4 : GACE on the binary image

	GLACE for Grouped Local Automated Cell Extractor
	

	# Paths d'entree 
	path_fused_0 : image fusionnee a l'instant t (celle dont on connait la segmentation)
	path_fused_1 : image fusionnee a l'instant t+1 (celle qu'on souhaite reconstruire localement)
	path_seg_0 : image segmentee a l'instant t
	path_membrane_prefix+'.[ext|theta|phi].inr' (optionel) : paths des images de maxima (ext) de membranes rehaussees et de leurs angles associes (theta, phi), afin d'eviter de les recalculer
	path_vector (optionel) : path vers le champ de deformation calcule (via blockmatching) entre t (flo) et t+1 (ref) : T_flo<-ref
							 Si le path existe deja, le champ de deformation n'est pas recalcule, sinon il est calcule et conserve

	# Label d'interet de l'instant 0, definissant la zone d'etude propagee a l'instant 1
	labels_of_interest : liste d'entiers correspondant aux labels de l'instant t dont on souhaite reconstruire le signal de membrane a t+1. Eventuellement, on peut ne donner qu'un entier (equivalent de LACE).
						Defaut : 'all' (signifie qu'on propage tous les labels qui existent dans l'image path_seg_0). Si 'all', le parametre background fixe le ou les label(s) de fond (defaut : [0,1]).


	# Path de sortie 
	path_membrane_prefix+'.[ext|theta|phi].inr' (optionel) : paths de sauvegarde des images de maxima (ext) de membranes rehaussees et de leurs angles associes (theta, phi)
	path_vector (optionel) : cf paths d'entree
	path_bin (optionel) : path de sauvegarde de l'image des membranes binarisees (image envoyee en entree de l'etape de vote de tenseurs)
	path_output (optionel) : path de sauvegarde de l'image reconstruite de sortie (par defaut : None)

	# Mask parameters
	rayon_dil=3.6 (defaut, exprime en reelles) : rayon de dilatation de la zone d'etude propagee a partir de l'instant t vers l'instant t+1

	# Membrane reconstruction parameters
	sigma_membrane=0.9 (defaut, unites reelles, adapte a des images de telles que Patrick/Ralph/Aquila) : parametre de rehaussement des 
                                                                                   membranes avant binarisation de celles-ci
	sensitivity=0.99 (defaut) : parametre de calcul des seuils anisotropiques selon un critere de sensibilite 
                                (true positive rate) : seuil = #(classe membrane>=seuil)/#(classe membrane) 
	manual=False (defaut) : parametre activant le mode manuel (ie parametrise) du seuil de binarisation des membranes si egal a True
	manual_sigma=7 (defaut) : parametre manuel d'initialisation pour le calcul du seuil de binarisation des membranes

	hard_thresholding=False (defaut) : Si echec de la precedente methode de seuillage des membranes par calcul automatique de seuils directionnels,
                              possibilite de choisir un seuillage dur global sur l'image en mettant cette option a True
	hard_threshold=1.0 (defaut) : Si hard_thresholding est a True, seuillage des membranes rehaussees via ce seuil 
	                          (1.0 : adaptee pour le time-point t001 d'Aquila par exemple)

	sigma_TV=3.6 (defaut, reelles, adapte a des images de resolution 0.3um^3) : parametre de propagation des membranes 
                                                                                par votes de tenseurs 
	sigma_LF=0.9 (defaut, reelles) : parametre de lissage gaussien de l'image des membranes reconstruite
	sample=0.2 (defaut) : echantillonne les votants de l'image initiale selon le coefficient (influe sur la vitesse de traitement)

	# Conserver ou non certaines images intermediaires
	keep_vector=False : conserver la transformation non lineaire T_t<-t+1 permettant de transformer l'image a t sur l'image a t+1 (se met a True automatiquement si path_vector est renseigne)
	keep_membrane=False : conserver toutes les images de l'etape de membranes (image d'extrema et d'angles) (se met a True automatiquement si path_membrane_prefix est renseigne)
	keep_all=False : conserver toutes les images intermediaires servant au calcul
	
	# Divers
	verbose=False : verbosite de la fonction, affiche par exemple les fichiers generes temporairement et ceux qui sont supprimes
	
	# Parallelisme 
	nb_proc=7 : nombre de processus lances en parallele pour LACE 

	'''

	# Step 1 : non-linear registration between path_fused_0 and path_fused_1

	# Step 2 : membrane reinforcement on global image

	# Step 3 : init and build the binary image (loop of short LACEs)

		# Step 3.1 : LACE on the label

		# Step 3.2 : adding to the main binary image

	# Step 4 : GACE on the binary image

	if type(labels_of_interest)==int:
		labels_of_interest=[labels_of_interest]



	# Test existence des images d'entree
	assert os.path.exists(path_fused_0), 'Miss file '+path_fused_0
	assert os.path.exists(path_fused_1), 'Miss file '+path_fused_1
	assert os.path.exists(path_seg_0), 'Miss file '+path_seg_0

	# Multi process import
	from multiprocessing import Process, Queue, Pool

	# Parametres de conservation des images generees
	if not keep_membrane:
		keep_membrane=keep_all
	if not keep_vector:
		keep_vector=keep_all

	# Definition des paths d'images intermediaires
	path_WORK=''
	if path_output:
		path_WORK = os.path.dirname(path_output).rstrip(os.path.sep)+os.path.sep
	else:
		path_WORK = os.path.dirname(path_fused_0).rstrip(os.path.sep)+os.path.sep
	if path_WORK==os.path.sep:
		path_WORK=os.path.curdir+os.path.sep

	tmp_ID = random_number()	

	keep_output=True
	if not path_output:
		keep_output=False
		path_output = path_WORK+'tmp_output_'+tmp_ID+'.inr'

	path_affine_trsf = path_WORK+'tmp_affine_'+tmp_ID+'.txt'

	if path_vector:
		keep_vector = True
	else:
		path_vector = path_WORK+'tmp_vectorfield_'+tmp_ID+'.inr'
	if path_membrane_prefix:
		keep_membrane = True
	else:
		path_membrane_prefix=path_WORK+'tmp_membrane_'+tmp_ID+''

	path_seg_trsf = path_WORK+'seg_trsf_at_1_'+tmp_ID+'.inr'
	path_tmp = path_WORK+'tmp_threshold_'+tmp_ID+'.inr'

	### Path de sortie ###
	if not os.path.isdir(path_WORK):
		try:
			os.mkdir(path_WORK)
		except Exception :
			print "Unexpected error: unable to create working directory"

	rayon_dil_voxel=None

	if rayon_dil:	# dilation of all the bounding boxes
		rayon_dil_voxel = rayon_dil / imread(path_fused_1).voxelsize[0]
		rayon_dil_voxel = int(rayon_dil_voxel+0.5)

	### Calculs ###

	# Step 1


	if not os.path.exists(path_vector):
		non_linear_registration(path_fused_0, path_fused_1, '/dev/null', path_affine_trsf, '/dev/null', path_vector, verbose=verbose)

	# Projection de la zone d'interet a l'instant t+1
	apply_trsf(path_seg_0, path_trsf=path_vector, path_output=path_seg_trsf, template=path_fused_1, nearest=True, verbose=verbose)


	# Step 2

	if (not os.path.exists(path_membrane_prefix+".ext.inr")) or (not os.path.exists(path_membrane_prefix+".theta.inr")) or (not os.path.exists(path_membrane_prefix+".phi.inr")):
		# Renforcement des membranes de l'image fusionnee a l'instant t+1 + extraction des maxima directionnels
		membrane_renforcement(path_fused_1, prefix_output=path_membrane_prefix, path_mask=None,  init=sigma_membrane, verbose=verbose)

	# Step 3 
	bboxes=boudingboxes(path_seg_trsf, verbose=verbose)
	if rayon_dil:	# dilation of all the bounding boxes
		rayon_dil_voxel = rayon_dil / imread(path_fused_1).voxelsize[0]
		rayon_dil_voxel = int(rayon_dil_voxel+0.5)
		dilation_tuple = (0, -rayon_dil_voxel, -rayon_dil_voxel, -rayon_dil_voxel, rayon_dil_voxel, rayon_dil_voxel, rayon_dil_voxel)
		import operator
		for x, b in bboxes.iteritems():
			b=map(operator.add, b, dilation_tuple) 
			if b[1] < 0:	# origine en x >= 0
				b[1] = 0
			if b[2] < 0:	# origine en y >= 0
				b[2] = 0
			if b[3] < 0:	# origine en z >= 0
				b[3] = 0
			# NB : les coordonnees du point "final" de la bounding box peuvent aller au dela de la dimension de l'image d'origine sans que cela n'affecte la suite du programme
			bboxes[x]=tuple(b)

	if type(labels_of_interest)==str and labels_of_interest=='all':
		#from lineage_test import labelsInImage
		#labels_of_interest=labelsInImage(imread(path_seg_0), background=background)
		labels_of_interest=[x for x in bboxes.keys() if not background.count(x)]

	pool=Pool(processes=nb_proc)
	mapping=[]

	for label_of_interest in labels_of_interest:

		# Step 3.1
		#path_local_bin=path_membrane_prefix+'.'+str(label_of_interest)+'.inr'

		# 	light_LACE METHOD INTERFACE (since ASTEC-170327):
		#		path_mask, label_of_interest, bbox, path_membrane_prefix, path_bin, rayon_dil, sigma_membrane, manual, manual_sigma, hard_thresholding, hard_threshold, sensitivity, verbose=parameters
		#parameters=(path_seg_trsf, label_of_interest, bboxes[label_of_interest], path_membrane_prefix, path_local_bin, \
		parameters=(path_seg_trsf, label_of_interest, bboxes[label_of_interest], path_membrane_prefix, None, \
			rayon_dil, manual, manual_sigma, hard_thresholding, hard_threshold, sensitivity, \
			verbose)
		if verbose:
			print 'Running light_LACE(' + str(parameters) + ') ...'

		mapping.append(parameters)


	outputs=pool.map(light_LACE, mapping)
	pool.close()
	pool.terminate()    

	path_union_of_local_bins=path_membrane_prefix+'.bin_union.inr'
	createImage(path_union_of_local_bins, path_fused_1, '-o 1', verbose=verbose)

	for path_local_bin in outputs:
		# Step 3.2
		label_of_interest = int(path_local_bin.split(os.path.sep)[-1].split('.')[-2])
		bbox=bboxes[label_of_interest]
		patchLogic(path_local_bin, path_union_of_local_bins, path_union_of_local_bins, bbox, Mode='or', verbose=verbose)
		#Logic(path_union_of_local_bins, path_local_bin, path_union_of_local_bins, Mode='or', verbose=verbose)
		if not keep_all:
			cmd='rm ' + path_local_bin
			print cmd
			os.system(cmd)

	# Suppression fichiers temporaires
	files_to_rm = path_affine_trsf + ' '
	if not keep_vector:
		files_to_rm += path_vector + ' '
	if files_to_rm:
		if verbose:
			print "Deleting temporary files: \n" + files_to_rm
		os.system("rm " + files_to_rm)

	keep_union_of_local_bins=False
	if path_bin:
		if path_bin != path_union_of_local_bins:
			keep_union_of_local_bins=False
			copy(path_union_of_local_bins, path_bin, verbose)
		else:
			keep_union_of_local_bins=True

	# Step 4

	GACE(path_union_of_local_bins, binary_input=True, path_output=path_output, 
		sigma_membrane=sigma_membrane, manual=manual, manual_sigma=manual_sigma, sensitivity=sensitivity, hard_thresholding=hard_thresholding, hard_threshold=hard_threshold, sigma_TV=sigma_TV, sigma_LF=sigma_LF, sample=sample,
		 keep_membrane=True, keep_all=False, verbose=verbose)

	reconstructed_image_1=imread(path_output)

	# Suppression fichiers temporaires
	files_to_rm=""

	if os.path.exists(path_membrane_prefix+'.hist.txt'):
		files_to_rm += path_membrane_prefix+'.hist.txt '

	if os.path.exists(path_membrane_prefix+'.bin.inr'):
		files_to_rm += path_membrane_prefix+'.bin.inr '

	if os.path.exists(path_seg_trsf):
		files_to_rm += path_seg_trsf + ' '

	if not keep_membrane:
		files_to_rm += path_membrane_prefix+'.ext.inr '
		files_to_rm += path_membrane_prefix+'.theta.inr '
		files_to_rm += path_membrane_prefix+'.phi.inr '
	if not keep_union_of_local_bins:
		files_to_rm += path_membrane_prefix+'.bin_union.inr '

	if not keep_output and os.path.exists(path_output):
		files_to_rm += path_output

	if files_to_rm:
		if verbose:
			print "Deleting temporary files: \n" + files_to_rm
		os.system("rm " + files_to_rm)


	return reconstructed_image_1
#
#	Calcul du flux recu par un cylindre
#		calcul integral pour le flux direct
#		modelisation anisotrope pour le flux diffus, d'apres la methode de (Perez et. al 1990)
#
#	El Mehdi Hamdani et Edouard Walther
#	AREP, 16 avenue d'Ivry
#
import math
global pi
pi = 3.141592653589793

# Cette fonction permet de retourner la valeur du flux diffus circumsolaire, a partir des valeurs du flux diffus horizontal, du flux direct normal, de la hauteur solaire et de l'angle d'incidence
def flux_diffus_circumsolaire_perez(diffus_horizontal, direct_normal, hauteur_solaire, angle_incidence):

	# Valeur de la constante solaire (flux solaire moyen qui parvient sur la terre a la limite de l'atmosphere) en W/m²
	I0 = 1368
	# angle solaire zenithal
	ksi_z = 90 - math.degrees(hauteur_solaire)
	# composante decrivant la transparence du ciel
	epsilon = (((diffus_horizontal + direct_normal)/diffus_horizontal) + 5.535 * 0.000001 * math.pow(ksi_z, 3))/(1 + 5.535 * 0.000001 * math.pow(ksi_z, 3))
	# Masse atmospherique (Kasten and Young 1989)
	a = max(0, math.cos(angle_incidence))
	b = max(math.cos(math.radians(85)), math.cos(math.radians(ksi_z)))
	ma = 1/(math.cos(math.radians(ksi_z)) + 0.50572 * ((96.07995 - ksi_z)**-1.6364))
	

	# Delta represente l'eclairement du ciel
	delta = ma * diffus_horizontal / I0
	
	# coefficients donnes par PEREZ et al (1990)
	if epsilon <= 1.065:
		f11 = 0.013
		f12 = 0.764
		f13 = -0.1
		f21 = -0.058
		f22 = 0.127
		f23 = -0.023
	elif epsilon <= 1.23:
		f11 = 0.095
		f12 = 0.920
		f13 = -0.152
		f21 = 0
		f22 = 0.051
		f23 = -0.02
	elif epsilon <= 1.5:
		f11 = 0.464
		f12 = 0.421
		f13 = -0.28
		f21 = 0.064
		f22 = -0.051
		f23 = -0.002
	elif epsilon <= 1.95:
		f11 = 0.759
		f12 = -0.009
		f13 = -0.373
		f21 = 0.201
		f22 = -0.382
		f23 = 0.01
	elif epsilon <= 2.8:
		f11 = 0.976
		f12 = -0.4
		f13 = -0.436
		f21 = 0.271
		f22 = -0.638
		f23 = 0.051
	elif epsilon <= 4.5:
		f11 = 1.176
		f12 = -1.254
		f13 = -0.462
		f21 = 0.295
		f22 = -0.975
		f23 = 0.129
	elif epsilon <= 6.2:
		f11 = 1.106
		f12 = -1.563
		f13 = -0.398
		f21 = 0.301
		f22 = -1.442
		f23 = 0.212
	else:
		f11 = 0.934
		f12 = -1.501
		f13 = -0.271
		f21 = 0.42
		f22 = -2.917
		f23 = 0.249
	
	# coefficient ponderant les densites du flux circumsolaire
	K1 = max(0, f11 + f12 * delta + pi * ksi_z * f13 / 180)
	
	# coefficient ponderant les densites du flux de l'horizon
	K2 = f21 + f22 * delta + pi * ksi_z * f23 / 180
	
	# flux diffus provenant du rayonnement circumsolaire W/m²
	diffus_circumsolaire = diffus_horizontal * K1 * a/b
	
	#retourne la valeur du flux diffus circumsolaire
	return diffus_circumsolaire

def flux_diffus_horizon_voute_perez(diffus_horizontal, direct_normal, hauteur_solaire, angle_incidence):
# Cette fonction permet de retourner la valeur du flux diffus circumsolaire, a partir des valeurs du flux diffus horizontal, du flux direct normal, de la hauteur solaire et de l'angle d'incidence

	# Valeur de la constante solaire (flux solaire moyen qui parvient sur la terre a la limite de l'atmosphere) en W/m²
	I0 = 1368
	# angle solaire zenithal
	ksi_z = 90 - math.degrees(hauteur_solaire)
	# composante decrivant la transparence du ciel
	epsilon = (((diffus_horizontal + direct_normal)/diffus_horizontal) + 5.535 * 0.000001 * math.pow(ksi_z, 3))/(1 + 5.535 * 0.000001 * math.pow(ksi_z, 3))
	
	# Masse atmospherique (Kasten and Young 1989)
	a = max(0, math.cos(angle_incidence))
	b = max(math.cos(math.radians(85)), math.cos(math.radians(ksi_z)))
	ma = 1/(math.cos(math.radians(ksi_z)) + 0.50572 * ((96.07995 - ksi_z)**-1.6364))
	
	
	# Delta represente l'eclairement du ciel
	delta = ma * diffus_horizontal / I0
	
	if epsilon <= 1.065:
		f11 = 0.013
		f12 = 0.764
		f13 = -0.1
		f21 = -0.058
		f22 = 0.127
		f23 = -0.023
		
	elif epsilon <= 1.23:
		f11 = 0.095
		f12 = 0.920
		f13 = -0.152
		f21 = 0
		f22 = 0.051
		f23 = -0.02
		
	elif epsilon <= 1.5:
		f11 = 0.464
		f12 = 0.421
		f13 = -0.28
		f21 = 0.064
		f22 = -0.051
		f23 = -0.002
		
	elif epsilon <= 1.95:
		f11 = 0.759
		f12 = -0.009
		f13 = -0.373
		f21 = 0.201
		f22 = -0.382
		f23 = 0.01
		
	elif epsilon <= 2.8:
		f11 = 0.976
		f12 = -0.4
		f13 = -0.436
		f21 = 0.271
		f22 = -0.638
		f23 = 0.051
		
	elif epsilon <= 4.5:
		f11 = 1.176
		f12 = -1.254
		f13 = -0.462
		f21 = 0.295
		f22 = -0.975
		f23 = 0.129
		
	elif epsilon <= 6.2:
		f11 = 1.106
		f12 = -1.563
		f13 = -0.398
		f21 = 0.301
		f22 = -1.442
		f23 = 0.212
		
	else:
		f11 = 0.934
		f12 = -1.501
		f13 = -0.271
		f21 = 0.42
		f22 = -2.917
		f23 = 0.249
		
	# coefficients donnes par PEREZ et al (1990)
	# coefficient ponderant les densites du flux circumsolaire
	K1 = max(0, f11 + f12 * delta + pi * ksi_z * f13 / 180)
	
	# coefficient ponderant les densites du flux de l'horizon
	K2 = f21 + f22 * delta + pi * ksi_z * f23 / 180
	
	# flux diffus provenant de l'horizon et de la voute celeste W/m²
	diffus_horizon_voute = diffus_horizontal * ((1 - K1) * (0.5 * (1 + math.cos(pi/2))) + K2 * math.sin(pi/2))
	
	# retourne la valeur du flux diffus de l'horizon et de la voute celeste
	return diffus_horizon_voute
	
# Cette fonction permet de retourner la valeur du flux solaire absorbe par un individu, a partir des valeurs des flux solaires diffus et direct horizontaux, de l'albedo du sol, du jour de l'annee, de l'heure, du numero du fuseau, de la latitude et de la longitude du site
def flux_waldhani (diffus_horizontal, direct_horizontal, albedo, jour, heure, fuseau, latitude, longitude_est): 

	# source : Wikipedia. L'equation du temps permet de calculer le temps solaire vrai.
	equation_temps = 7.678 * math.sin(1.374 + (2 * pi * (jour - 81)/365)) - 9.87 * math.sin(2 * (2 * pi * (jour - 81)/365)) 
	
	# le temps solaire vrai permet de calculer la position du soleil
	temps_solaire_vrai = heure + (longitude_est/15) - fuseau - (equation_temps/60) 
	
	# l'unite de l'angle horaire doit être en radians car dans le langage python l'angle des fonctions trigonometriques est considere comme etant en radians
	angle_horaire = math.radians((temps_solaire_vrai - 12) * 15)
	
	# declinaison en radians
	declinaison = math.asin(0.4 * math.sin(math.radians(0.986 * jour - 80)))
	
	# hauteur solaire en radians
	hauteur_solaire = math.asin(math.sin(math.radians(latitude)) * math.sin(declinaison) + math.cos(math.radians(latitude)) * math.cos(declinaison) * math.cos(angle_horaire))
	
	#Afin d'eviter des erreurs dans la suite du calcul nous limitons certaines valeurs de hauteur solaire a une valeur nulle
	if hauteur_solaire < math.radians(2):
		hauteur_solaire = 0
	
	#L'azimut est l'angle compte a partir du sud, positivement vers l'ouest, negativement vers l'est (temps_solaire_vrai < 12), entre le plan vertical du soleil a l'instant donne et le plan meridien local.
	if temps_solaire_vrai < 12:
		azimut = -1 * math.acos((math.sin(math.radians(latitude)) * math.cos(declinaison) * math.cos(angle_horaire) - math.cos(math.radians(latitude)) * math.sin(declinaison)) / math.cos(hauteur_solaire)) 
	else:
		azimut = math.acos((math.sin(math.radians(latitude)) * math.cos(declinaison) * math.cos(angle_horaire) - math.cos(math.radians(latitude)) * math.sin(declinaison)) / math.cos(hauteur_solaire))
	
	
	#Cette formule permet de calculer le cosinus de l'angle d'incidence
	#	(angle entre la normale au plan recepteur et le rayon solaire incident)
	# Le math.cos(0) signifie que le plan suit la trajectoire du soleil.
	cosinus_angle_incidence = math.cos(hauteur_solaire)*math.sin(0.5*pi)*math.cos(0)+math.sin(hauteur_solaire)*math.cos(0.5*pi)

	#pour eviter une division par 0 dans la formule ligne 224, lorsque la hauteur solaire est egale a 0 on considere que le flux direct normal est nul	
	if hauteur_solaire == 0:
		direct_normal = 0
	else:
		direct_normal = direct_horizontal/math.sin(hauteur_solaire)

	if diffus_horizontal == 0:
		flux_diffus_circumsolaire = 0
		flux_diffus_voute_horizon = 0	
	else:
		flux_diffus_circumsolaire = flux_diffus_circumsolaire_perez(diffus_horizontal, direct_normal, hauteur_solaire, math.acos(cosinus_angle_incidence))
		# calcul de l'azimut de la face arriere de la paroi verticale en degre
		if math.degrees(azimut) < 0:
			azimut_arr = math.degrees(azimut) + 180
		else:
			azimut_arr = math.degrees(azimut) - 180
		
		azimut_arriere = math.radians(azimut_arr)
		#cosinus de l'angle d'incidence de la face arriere de la paroi verticale
		cosinus_angle_incidence_arriere = math.cos(hauteur_solaire) * math.sin(0.5 * pi) * math.cos(azimut - azimut_arriere) + math.sin(hauteur_solaire) * math.cos(0.5 * pi)
		
		flux_diffus_voute_horizon = flux_diffus_horizon_voute_perez(diffus_horizontal, direct_normal, hauteur_solaire, math.acos(cosinus_angle_incidence)) + flux_diffus_horizon_voute_perez(diffus_horizontal, direct_normal, hauteur_solaire, math.acos(cosinus_angle_incidence_arriere))
	
	if flux_diffus_circumsolaire < 0:
		flux_diffus_circumsolaire = 0

	if flux_diffus_voute_horizon < 0:
		flux_diffus_voute_horizon = 2 * diffus_horizontal * 0.5
	
	if hauteur_solaire == 0:
		flux_direct = 0
	else:
		flux_direct = direct_horizontal * cosinus_angle_incidence / math.sin(hauteur_solaire)
		
	flux_reflechi_direct = albedo * direct_horizontal * cosinus_angle_incidence	* (1 - math.cos(pi/2))/2
	flux_reflechi_diffus = albedo * diffus_horizontal * (1 - math.cos(pi/2))/2
	
	# geometrie du cylindre
	rayon_cylindre = 0.17
	hauteur_cylindre = 1.73
	# calcul des flux recus par le cylindre 
	# 	> direct = anisotrope (calcul integral)
	flux_direct_cylindre = 2 * flux_direct * rayon_cylindre * hauteur_cylindre
	# 	> diffus = anisotrope circumsolaire (calcul integral) + isotrope diffus et horizon
	flux_diffus_cylindre = (2 * flux_diffus_circumsolaire + pi * flux_diffus_voute_horizon) * rayon_cylindre * hauteur_cylindre		
	# 	> reflechi  = anisotrope (calcul integral) + isotrope
	flux_reflechi_cylindre = (2 * flux_reflechi_direct + 2 * pi * flux_reflechi_diffus) * rayon_cylindre * hauteur_cylindre	
		
	#Valeurs en Watts
	return(flux_direct_cylindre, flux_diffus_cylindre, flux_reflechi_cylindre)

def temperature_mrt_fluxsolaire(temperature_mrt_classique, diffus_horizontal, direct_horizontal, albedo, jour, heure, fuseau, latitude, longitude_est):
	sigma = 5.67 * 0.00000001 #constante de stefan boltzmann
	epsilon = 0.97#emissivite du corps humain
	#temperature radiante en Kelvin
	temperature_radiante = math.pow(temperature_mrt_classique + 273.15, 4)
	f_eff = 0.75 #facteur de surface de rayonnement effectif pour un individu debout d'après JENDRITZKY and NUBLER
	alpha_clo = 0.8 #coefficient d'absorption des vêtements
	# geometrie cylindre
	rayon_cylindre = 0.17
	hauteur_cylindre = 1.73
	S_cylindre = pi * rayon_cylindre * hauteur_cylindre
	
	#flux direct, diffus, reflechi Waldhani en W/m²
	flux_direct, flux_diffus, flux_reflechi = flux_waldhani(diffus_horizontal, direct_horizontal, albedo, jour, heure, fuseau, latitude, longitude_est) / S_cylindre
	temperature_mrt_solaire = math.pow(temperature_radiante + (alpha_clo * flux_direct + f_eff * alpha_clo * flux_diffus + f_eff * alpha_clo * flux_reflechi)/(sigma * epsilon), 0.25) - 273.15
	
	return temperature_mrt_solaire
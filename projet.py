import sys   #sys pour pouvoir recuperer notre commande à la fin du code

#Q1 comptez le nombre de reads mappés : 

def countMappedReads(sam_file):  
    mapped_count = 0  # on initialise notre compteur à 0

    with open(sam_file, 'r') as file: #on ouvre le fichier sam_file et on le place ds le parametre file/ with as c pour ouvrire le fichier ds la variable file et le 'r' sert à lire le fichier si on avait mit w on pourrait l'ecrire
        for line in file:
            if line.startswith('@'):  #on ignore toute ligne qui commence par @ (l'entete)
                continue 
            
            columns = line.strip().split('\t') 
            #on eneleve les espaces du debut et de la fin avec stript  de la chaine de caractere 
            #la colonne devient une liste de chaine de caractere  couper à chaque tab puisque SAM est un fichier tabulé
            flag = int(columns[1])  #on transforme la chaine de caractères de la colonne flag en un entier pour pouvoir traiter les bits
            
            if (flag & 4) == 0: #read est mappé si 3ème bit de droite est à 0/flag=4 est non mappé puisque en bineaire le bit3 est 1 
                mapped_count += 1  #on rajoute 1 au compteur

    print(f"Nombre total de reads mappés : {mapped_count}") #f pour qu'il interprete la valeur à l'interieur des {}
    

#Q2 comptez le nombre de raeds par flag : faut traduire le flag ici et resumer uniquement ceux qui ont un sens 
def countFlags(sam_file):
    # on crée un dictionnaire pour compter les cas où les deux paires de reads sont mappés ou un seul des deux est mappé
    paired_status = {
        "both_mapped": 0,   # on crée une clé pour les deux reads de la paire qui sont mappés et on initialise la vlaur de la clé à 0
        "one_mapped": 0,    # on crée une clé pour un seul read de la paire est mappé et on initialise la vlaur de la clé à 0
    }

    read_pairs = {}  # on crée un dictionnaire vide pour stocker les flags des reads par leur QNAME car une paire de reads ont le meme Qname

    with open(sam_file, 'r') as file: #on ouvre le fichier sam_file et on le place ds le parametre file/ with as c pour ouvrire le fichier ds la variable file et le 'r' sert à lire le fichier si on avait mit w on pourrait l'ecrire
        for line in file:
            if line.startswith('@'): 
                continue  #ici il lis ligne par ligne du c=fichier sam et à chaque fois u'il trouve un @ il ignore la ligne car c'est les lignes d'en tete qui commencent par un @
            
            columns = line.strip().split('\t')  #on eneleve les espaces du debut et de la fin avec strip  de la chaine de caractere et la colonne devient une liste de chaine de caractere  couper à chaque tab puisque SAM est un fichier tabulé
            qname = columns[0]  
            flag = int(columns[1])  # on transforme le FLAG en entier

            
            if qname not in read_pairs:  # ici on va regrouper les flags des reads ayant le même QNAME (foward et reverse) et on crée une condition pour vérifier si le qname existe déjà comme clé dans le dictionnaire read_pairs/ S'il n'existe pas ca va etre la première fois qu'on rencontre ce read.
                read_pairs[qname] = [] # dans le dictionnaire read_pairs on crée une liste vide associé à la clé qui est le qname rencontré
            read_pairs[qname].append(flag) # grace à append on peut ajouter le flag à la fin de la liste qui porte comme clé le qname rencontré

    # Analyser les paires de reads
    for qname, flags in read_pairs.items(): # pour chaque qname on verifie si les reads sont paires (si le reverse et le foward existe) et on ignorer les cas où il n'y a pas deux reads pour une paire
        if len(flags) != 2: #si le nombre de flag pour chaque qname est pas egal à 2 on ignore ce qname 
            continue
        
        # on vérifie si les deux reads de la paire qui portent le meme Qname sont mappés en créant deux variables et en faisant une opération booleenne entre le flag du read et le bit 4        
        first_mapped = (flags[0] & 0x4)==0  #  le premier read est mappé si le flag du premier read & 0x4 == 0 
        second_mapped = (flags[1] & 0x4)==0  # Le second read est mappé si le flag du second read & 0x4 == 0

        if first_mapped and second_mapped: #ici on va faire une condition pour savoir si les deux reads sont mappées et rajouter un curseur compter +1 quand c'est le cas
            paired_status["both_mapped"] += 1 #si la condistion booleene renvoie vraie pour les deux reads mappées on rajouté plus 1 dans le curseur de la clé du dictionnaire paired statut
        elif first_mapped or second_mapped: # si non si la condition booleenne renvoie vraie uniquement un seul des read est mappé on rajoute un curseur à la clé one mapped du dictionnaire paired statut 
            paired_status["one_mapped"] += 1

    # Afficher les résultats
    print("Comment les reads et paires de reads :")
    print(f"Paires avec les deux reads mappés : {paired_status['both_mapped']}") # f c'est pour qu'il prenne la valeur qu'il y a entre les {} et indique qu'il doit afficher le contenu des clés du dictionnaire en format string(chaines de caractère)
    print(f"Paires avec un seul read mappé : {paired_status['one_mapped']}")

#Q3: Où les reads sont-ils mappés ? L'alignement est-il homogène le long de la séquence de référence ? 
def countReadsBySegment(sam_file, segment_size=10000): #On crée une fonction pour compter les reads dont la position de début se situe sur un segement du chromosome de taille predéfinie (ici par defaut 10000pb) . cette fonction prend en entrée le fichier SAM et une taille de segment par défaut de 10000 bases
    chromosome_segments = {} # ici on crée un dicitonaire vide pour stocker les segments de taille de 10000 pb du chromosome et le nombre de reads associées à chaque segment 


    with open(sam_file, 'r') as file: # on ouvre le fichier sam en mode lecture grace au 'r'
        for line in file: #pour chaque ligne du fichier qu'il va analyser une par une 
            if line.startswith('@'):  #on crée une condition si la ligne commence par un @
                continue #il l'ignore on va donc ignorer les lignes d'en tete qui commencent par un @

            columns = line.strip().split('\t') # ensuite on supprime les espace qu'on peut avoir au debut et à la fin des lignes avec strip et oavec split on sépare la ligne en colonnes en utilisant le caractère de tabulation comme séparateur
            chromosome = columns[2]  # on designe que le nom du chromosome est dans la 3ème colonne (index 2)
            start_position = int(columns[3])  # ici on désigne que la position de debut du read se trouve dans la colonne 2 ayant l'indice 3 et on transforme le chiffre en entier avec 'int' car dnas le fichier c'est une chaine de caractère

            segment_index = start_position // segment_size # Calcul l'indice du segment auquel appartient le read en divisant la position de départ par la taille du segment

        
            if chromosome not in chromosome_segments: # on fait une condition pour vérifier si le chromosome n'existe  aps deja dans le dictionnaire 
                chromosome_segments[chromosome] = {}  # si non on ajoute une entrée pour ce chromosome avec un dictionnaire vide pour ses segments           
            if segment_index not in chromosome_segments[chromosome]: #on crée une deuxième condition pour verifier si l'index du segment existe, s'il n'existe pas on initialise le compteur dans la prochaine ligne 
                chromosome_segments[chromosome][segment_index] = 0  # on initialise le compteur de reads pour ce segment à 0

            chromosome_segments[chromosome][segment_index] += 1  # le compteur ajoute +1

    print("Distribution des reads par segment :") # on affiche le resultats de la distibution des reads par segment de 10000 pb 
    print("Chromosome\tSegment_Start\tReads_Count") # on affiche une en-tête avec les colonnes : Chromosome, Début du segment, Nombre de reads
    for chromosome, segments in chromosome_segments.items(): #pour chaque chromosome dans le dictionnaire chromosome_segments
        for segment_index, count in sorted(segments.items()): #on trie les segments par index avec la fonction 'sorted' pour les afficher dans l'ordre
            segment_start = segment_index * segment_size     #ici il va compter la position de debut du segement en multipliant l'index du segment x10000pb
            print(f"{chromosome}\t{segment_start}\t{count}") # on affiche les resultats pour chaque chromosome la position de debut et le nombre de reads associé à ce segment


# Q4: Comptez le nombre de reads par tranche de scores mapQ : ici il faut faire un diagramme pour voir ce que ca represente et les faire par tranche de 10
def countReadsByMapQ(sam_file):  # on crée une fonction qui prend en entrée le fichier SAM et compte le nombre de reads pour chaque score MAPQ.
    mapq_counts = {} #on cree un dictionnaire vide

    with open(sam_file, 'r') as file: #on ouvre le fichier sam et le lire grace au 'r' comme un fichier 
        for line in file: # on crée  une boucle 'pour' où pour chaque ligne du fichier 
            if line.startswith('@'): #on crée une condition si la ligne commence par un @ on l'ignore (les lignes d'en tete)
                continue

            columns = line.strip().split('\t') # ensuite on supprime les espace qu'on peut avoir au debut et à la fin des lignes avec strip et oavec split on sépare la ligne en colonnes en utilisant le caractère de tabulation comme séparateur
            mapq = int(columns[4])  # on designe que la qualité de mapping est dans la 5ème colonne (index 4)

            if mapq not in mapq_counts: #on cree une condition si le mapq rencontré n'est pas dans le dictionnaire map_counts
                mapq_counts[mapq] = 0 # On initialise le compteur pour ce score MAPQ à 0

            mapq_counts[mapq] += 1 # on ajoute +1 au compteur

    print("Distribution des scores de qualité de mapping (MAPQ) :") #on affiche les résultats de la distribution des scores de qualité de mapping 
    for mapq, count in sorted(mapq_counts.items()):# grace à 'sorted' on trie les scores MAPQ dans l'ordre croissant
        print(f"Score MAPQ {mapq}: {count} reads") # on affiche le score et le nombre de read pour chaque score et 'f' prends la valeur compté entre les {}





# python3 projet.py mapping.sam   ->   sys.argv[0] = projet.py   sys.argv[1] = mapping.sam
        
input_sam = sys.argv[1]  # on récupérère le nom du fichier SAM en argument de la commande
countMappedReads(input_sam) # ici on appelle la fonction countMappedReads en lui passant le fichier sam (input_sam) comme argument.
countFlags(input_sam) # ici on appelle la fonction countFlags en lui passant le fichier sam (input_sam) comme argument.
countReadsBySegment(input_sam) # ici on appelle la fonction coutReadsBySegments en lui passant le fichier sam (input_sam) comme argument.
countReadsByMapQ(input_sam) # ici on appelle la fonction countReadsByMapq en lui passant le fichier sam (input_sam) comme argument.

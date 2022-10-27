##### Lancer une tache sur un cluster
1. Connexion ssh
``` 
ssh dahu
```
2. Création du script de soumission (voir exemple)
``` 
vim script_name.sh
```
3. Rendre le script exécutable
``` 
chmod +x script_name.sh
```
4. Soumission de la tâche
``` 
oarsub -S ./script_name.sh
```
Avant de lancer une tache sur le cluster pour la première fois, il faut s'assurer qu'il est possible de cloner le répertoire git sur le serveur. 

##### Récupérer des fichiers présents sur le cluster
Depuis sa machine
``` 
scp dahu:file_location_on_dahu/file_name file_location_on_local
```
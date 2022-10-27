#### Lancer une tache sur un cluster
Exécuter depuis la frontale du cluster :
``` 
oarsub -S script_name.sh
```

#### Récupérer des fichiers présents sur le cluster
Exécuter depuis sa machine (par exemple)
``` 
scp dahu:file_location_on_dahu/file_name file_location_on_local
```
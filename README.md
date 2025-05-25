TASK:
Реализовать один из простейших методов кластеризации (например, K-means) в виде
библиотечной функции. Подготовить датасет, на котором будет продемонстрирована
работоспособность алгоритма.

TO RUN:
gcc -o result_of_clustering main.c project_final.c
./result_of_clustering 
or
gcc -fPIC -c project_final.c
gcc -shared -o project_final.so project_final.o
gcc -o result_of_clustering main.c ./project_final.so

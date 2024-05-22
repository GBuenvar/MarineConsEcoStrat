Inspirado por las ideas de Pietronero Cristelli et al., M.A.Muñoz y Virginia Dominguez_García presentan el algoritom MusRank.

Este algoritmo sirve para la ordenacion de redes mutualistas de la manera mas nested posiible.
En Redes mutualistas, en las que se relacionan dos conjuntos de nodos, como pueden ser bienes exportados por paises, semillas
dispersadas por pajaros, etc., se pueden encontrar estructuras nested en las que los nodos especialistas de un conjuntos
se relacionan con los nodos generalistas del otro conjunto. En tal caso, definen uno de los conjuntos como el 
"activo" (países, dospersadores, peces) y otro como el "pasivo" (bienes, semillas, anemonas). El rol de cada conjunto
no es una mera cuestion de nomenclatura, sino que caracteriza la estructura de la red. 
Por ejemplo, en el caso de la red de exportaciones de productos por paises, se observa una ordenacion en la que los productos
que son exportados por pocos paises, son exportados por los paises que exportan muchos productos, mientras que los paises que 
exportan pocos productos exportan productos que son exportados por muchos paises. No hay especializacion en la red

- los paises generalistas exportan todo tipo de productos
- los productos generalistas son exportados por todo tipo de paises
- los productos que son exportados por pocos paises son exportados por los paises generalistas
- los paises que exportan pocos productos exportan productos generalistas

Este tipo de diseño confiere robusted a la red frente a perdidas de nodos, ya que los nodos generalistas pueden ser sustituidos
por otros nodos generalistas, mientras que los nodos especialistas no pueden ser sustituidos por otros nodos especialistas.


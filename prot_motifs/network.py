'''
NAME
    network.py

VERSION
    [1.0]

AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>

DESCRIPTION
    This module is part of the bio_analysis package.
    It works with regulatory networks to obtain some of its general traits such as connectivity and clustering

CATEGORY
    Regulatory networks

USAGE
    import bio_analysis.network

CLASSES
    net([network path]): defines a network as an object

FUNCTIONS
    num_nodes([name of network]): determines the number of nodes of the network
    num_edges([name of network]): determines the number of edges of the network
    get_edges([name of node]): obtains the edges of a certain node
    neighbors([name of node]): obtains the nodes that interact with a certain node
    grades([name of node]): obtains the node's grade
    clustering([name of node]): obtains the clustering coefficient of a certain node
    P_k([grade]): determines the number of nodes with a certain grade
    C_k([grade]): determines the average clustering of the nodes with a certain garde

GITHUB LINK
    https://github.com/daianna21/python_class/blob/master/bio_analysis/network.py

'''

# The class net is created
class net:
    # A particular object is instantiated with the path to the network
    def __init__(self, netpath):
       # Object attributes are defined: lists of nodes (genes), edges (interactions between genes) and the grades of
       # each node
        self.nodes = []
        self.edges= []

       #The file with the network is open: each line has the two interacting genes linked by an edge
        with open(netpath) as file:
            for line in file:
                genes=line.split()
                # Nodes are added in the list self.nodes
                if genes[0] not in self.nodes:
                    self.nodes.append(genes[0])
                if genes[1] not in self.nodes:
                    self.nodes.append(genes[1])
                #Edges are added in the list self.edged as tuples that contain the interacting genes names
                if (genes[0],genes[1]) not in self.edges:
                    self.edges.append((genes[0], genes[1]))
    #Methods of the object:
    def num_nodes(self):
        """
        This function gets the number of nodes of the network
                Parameters:
                        The object itself: the network
                Returns:
                        len(self.nodes) (int): The length of the list of unique nodes
        """
        #The list contains all genes only once
        return(len(self.nodes))

    def num_edges(self):
        """
        This function gets the number of edges of the network
                Parameters:
                       The object itself: the network
                Returns:
                       len(self.edges) (int): The length of the list of unique edges (tuples of genes)
        """
        #The list contains all existent interactions only once
        return(len(self.edges))

    def get_edges(self, node):
        """
        This function gets the interactions (edges) of a certain gene (node)
                Parameters:
                        The object itself: the network
                        node (str): the gene's name
                Returns:
                        edges (list): a list of all tuples that contain the gene
        """
        #Get all tuples that have gene's name
        edges=[tuple for tuple in self.edges if node in tuple]
        return(edges)

    def neighbors(self, node): #obtener vecinos de un nodo
        """
         This function gets the neighbors of a certain node: the genes with which a certain gene interacts
                 Parameters:
                         The object itself: the network
                         node (str): the gene's name
                 Returns:
                         neigh (list): a list of the interacting genes names
         """
        #Get all genes of the tuples in which the gene of interest is found
        neigh=[str for tuple in self.get_edges(node) for str in tuple if not str == node]
        return(neigh)

    def grades(self, node):
        """
         This function gets the grade of a certain node
                 Parameters:
                         The object itself: the network
                         node (str): the gene's name
                 Returns:
                         grade (int): the number of neighbors the node has
         """
        #The grade is given by the number of elements the list neigh of the gene has
        grade=len(self.neighbors(node))
        return(grade)

    def clustering(self, node):
        """
         This function calculates the clustering coefficient of a certain node
                 Parameters:
                         The object itself: the network
                         node (str): the gene's name
                 Returns:
                         clust (float): tha value of the clustering coefficient of a gene
         """
        # Clustering coefficient is given by the number of existent interactions between the neighbors of a node divided
        # by the total posible interactions between them:
        # Get the neighbors of the gene's neighbors that are neighbors of the gene as well, ignore self-regulation
        neig_neighbors= [tuple for tuple in self.edges if tuple[1] in self.neighbors(node)
                         and tuple[0] in self.neighbors(node) if not tuple[0] == tuple[1]]
        #Number of existent interactions
        len_neig_neighbors= len(neig_neighbors)
        #Number of total possible interactions: is divided by 2 since each pair of nodes is counted twice
        total_interact= (self.grades(node) * (self.grades(node) - 1)) / 2
        #Clustering is the ratio of the above values
        clust=len_neig_neighbors/total_interact
        return(clust)

    def P_k(self, grade):
        """
          This function calculates the number of nodes in the network with a certain grade (P(k))
                  Parameters:
                          The object itself: the network
                          grade (int): the grade to which connectivity wants to be determined
                  Returns:
                          connectivity (int): the number of genes that have that grade
        """
        #Get a list with all the genes that have that grade
        nodes_grade=[str for str in self.nodes if self.grades(str) == grade]
        #Get the number of genes the list has
        connectivity= len(nodes_grade)
        return(connectivity)

    def C_k(self, grade):
        """
          This function calculates the average clustering coefficient of all genes that have a certain grade (C(k))
                  Parameters:
                          The object itself: the network
                          grade (int): the grade to which average clustering wants to be determined
                  Returns:
                          aver_clust (float): the average clustering of the genes
        """
        #Get a list of the genes that have that certain grade
        nodes_grade=[str for str in self.nodes if self.grades(str) == grade]
        sum_clust=0
        #The clustering is calculated for each node of the list and added to sum_clust
        for nodo in nodes_grade:
          sum_clust=sum_clust + self.clustering(nodo)
        #The average is calculated
        aver_clust=sum_clust/len(nodes_grade)
        return(aver_clust)


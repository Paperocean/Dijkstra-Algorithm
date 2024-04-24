#include <iostream>
#include <limits>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

struct Target {
    int vertex;
    int distance;

    Target(int _vertex, int _distance) : vertex(_vertex), distance(_distance) {}
};

class ChangeablePriorityQueue {
private:
    vector<Target*> Q;
    int currentSize;

    int parent(int i) { return (i - 1) / 2; }
    int leftChild(int i) { return (2 * i + 1); }
    int rightChild(int i) { return (2 * i + 2); }

    void heapifyUp(int index) {
        while (index > 0 && Q[index]->distance > Q[parent(index)]->distance) {
            swap(Q[index], Q[parent(index)]);
            index = parent(index);
        }
    }

    void heapifyDown(int index) {
        int maxIndex = index;
        int left = leftChild(index);
        int right = rightChild(index);

        if (left < currentSize && Q[left]->distance > Q[maxIndex]->distance) {
            maxIndex = left;
        }

        if (right < currentSize && Q[right]->distance > Q[maxIndex]->distance) {
            maxIndex = right;
        }

        if (maxIndex != index) {
            swap(Q[index], Q[maxIndex]);
            heapifyDown(maxIndex);
        }
    }

public:
    ChangeablePriorityQueue() : currentSize(0) {}

    bool isEmpty() const {
        return currentSize == 0;
    }

    bool find(int temp) const {
        for (size_t i = 0; i < currentSize; i++) {
            if (Q[i]->vertex == temp)
                return true;
        }
        return false;
    }

    Target* peek() const {
        if (currentSize > 0)
            return Q[0];
        return nullptr;
    }

    void build(Target* target) {
        Q.push_back(target);
        heapifyUp(currentSize++);
    }

    void deleteMin() {
        if (isEmpty())
            return;

        Q[0] = Q[currentSize - 1];
        Q.pop_back();
        currentSize--;
        heapifyDown(0);
    }

    void decreaseKey(int id, int k) {
        if (id < 0 || id >= currentSize)
            return;

        if (k >= Q[id]->distance)
            return;

        Q[id]->distance = k;
        heapifyUp(id);
    }
};

struct AdjacencyMatrix {
    int** adjacencyMatrix;
    string* vertexLabels;
    int V;

    AdjacencyMatrix(int numVertices) : V(numVertices) {
        adjacencyMatrix = new int* [V];
        for (int i = 0; i < V; ++i) {
            adjacencyMatrix[i] = new int[V];
            for (int j = 0; j < V; ++j) {
                adjacencyMatrix[i][j] = -1;
            }
        }
        vertexLabels = new string[V];
    }

    ~AdjacencyMatrix() {
        for (int i = 0; i < V; ++i) {
            delete[] adjacencyMatrix[i];
        }
        delete[] adjacencyMatrix;
    }

    void setVertexLabel(int vertex, string label) {
		vertexLabels[vertex] = label;
	}

    string* getVertexLabels() {
        string* labels = new string[V];
        for (int i = 0; i < V; i++) {
			labels[i] = vertexLabels[i];
		}
        return labels;
    }

    void addEdge(int u, int v, int w) {
        adjacencyMatrix[u][v] = w;
    }

    void displayEdges(int vertex) {
        for (int i = 0; i < V; i++) {
            if (!vertexLabels[vertex].empty() && adjacencyMatrix[vertex][i] != -1) {
                cout << "Edge: " << vertexLabels[vertex] << " - " << vertexLabels[i] << " Weight: " << adjacencyMatrix[vertex][i] << endl;
            }
            else if (vertexLabels[vertex].empty() && adjacencyMatrix[vertex][i] != -1) {
				cout << "Edge: " << vertex << " - " << i << " Weight: " << adjacencyMatrix[vertex][i] << endl;
			}
		}
	}

    void display() {
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                if (adjacencyMatrix[i][j] == -1) {
                    cout << "INF ";
                }
                else {
                    cout << adjacencyMatrix[i][j] << " ";
                }
            }
            cout << endl;
        }
        cout << endl;
    }

    double dijkstra(int source, ChangeablePriorityQueue& pq) {
        auto start = chrono::high_resolution_clock::now();
        if (source < 0 || source >= V) {
            cout << "Invalid source vertex index!" << endl;
            return 0.0;
        }

        int* distance = new int[V];
        bool* visited = new bool[V];

        for (int i = 0; i < V; ++i) {
            distance[i] = numeric_limits<int>::max();
            visited[i] = false;
        }
        distance[source] = 0;

        pq.build(new Target(source, 0));

        // Algorytm Dijkstry
        while (!pq.isEmpty()) {
            Target* current = pq.peek(); // Wybieramy wierzchołek o najmniejszej odległości
            int u = current->vertex;
            pq.deleteMin(); // Usuwamy wybrany wierzchołek z kolejki

            if (visited[u]) continue; // Jeśli wierzchołek był już odwiedzony, przechodzimy do następnego

            visited[u] = true; // Oznaczamy wierzchołek jako odwiedzony

            // Aktualizujemy odległości dla sąsiadów wybranego wierzchołka
            for (int v = 0; v < V; ++v) {
                if (!visited[v] && adjacencyMatrix[u][v] != -1 && distance[u] != numeric_limits<int>::max() &&
                    distance[u] + adjacencyMatrix[u][v] < distance[v]) {
                    distance[v] = distance[u] + adjacencyMatrix[u][v];

                    int index = -1;
                    if (pq.find(v)) {
                        index = v;
                    }

                    if (index != -1) {
                        pq.decreaseKey(index, distance[v]); // Aktualizujemy odległość dla v
                    }
                    else {
                        pq.build(new Target(v, distance[v])); // Dodajemy v do kolejki | DLACZEGO DODAJEMY NA KONCU ??? PO CO??
                    } // Dodajemy lub aktualizujemy wierzchołek w kolejce
                }
            }
            delete current;
        }

        auto end = chrono::high_resolution_clock::now();

        // Wypisujemy obliczone odległości
       /* cout << "Distances from source vertex " << source << endl;
        for (int i = 0; i < V; ++i) {
            if (distance[i] == numeric_limits<int>::max()) {
                cout << "Vertex " << i << ": Unreachable" << endl;
            }
            else {
                cout << "Vertex " << i << ": " << distance[i] << endl;
            }
        }*/

        delete[] distance;
        delete[] visited; 
        return chrono::duration<double, milli>(end - start).count();
    }
};

struct AdjacencyList {
    pair<int, int>** adjacencyList;
    string* vertexLabels;
    int V;

    AdjacencyList(int numVertices) : V(numVertices) {
        adjacencyList = new pair<int, int>* [V];
        for (int i = 0; i < V; i++) {
            adjacencyList[i] = new pair<int, int>[V];
            for (int j = 0; j < V; j++) {
                adjacencyList[i][j] = make_pair(-1, -1);
            }
        }
        vertexLabels = new string[V];
    }

    ~AdjacencyList() {
        for (int i = 0; i < V; i++) {
            delete[] adjacencyList[i];
        }
        delete[] adjacencyList;
    }

    void setVertexLabel(int vertex, string label) {
        vertexLabels[vertex] = label;
    }

    string* getVertexLabel() {
		string* labels = new string[V];
        for (int i = 0; i < V; i++) {
			labels[i] = vertexLabels[i];
		}
		return labels;
	}

    void addEdge(int u, int v, int w) {
        adjacencyList[u][v] = make_pair(v, w);
        adjacencyList[v][u] = make_pair(u, w);
    }
    
    void displayEdges(int vertex) {
        for (int i = 0; i < V; i++) {
            if (!vertexLabels[vertex].empty() && adjacencyList[vertex][i].first != -1) {
                cout << "Edge: " << vertexLabels[vertex] << " - " << vertexLabels[i] << " Weight: " << adjacencyList[vertex][i].second << endl;
            }
            else if (vertexLabels[vertex].empty() && adjacencyList[vertex][i].first != -1) {
                cout << "Edge: " << vertex << " - " << adjacencyList[vertex][i].first << " Weight: " << adjacencyList[vertex][i].second << endl;
            }
		}
	}

    void display() {
        for (int i = 0; i < V; i++) {
            cout << "Vertex " << i << ": ";
            for (int j = 0; j < V; j++) {
                const auto& neighbor = adjacencyList[i][j];
                if (neighbor.first != -1 && neighbor.second == -1) {
                    cout << "(" << neighbor.first << ", INF) ";
                }
                else if (neighbor.first == -1) {
					cout << "INF ";
                }
                else {
                    cout << "(" << neighbor.first << ", " << neighbor.second << ") ";
                }
            }
            cout << endl;
        }
        cout << endl;
    }

    double dijkstra(int source, ChangeablePriorityQueue& pq) {
        auto start = chrono::high_resolution_clock::now();

        if (source < 0 || source >= V) {
            cout << "Invalid source vertex index!" << endl;
            return 0.0;
        }

        int* distance = new int[V];
        bool* visited = new bool[V];

        for (int i = 0; i < V; ++i) {
            distance[i] = numeric_limits<int>::max();
            visited[i] = false;
        }
        distance[source] = 0; 
        pq.build(new Target(source, 0)); 
        while (!pq.isEmpty()) {
            Target* current = pq.peek();
            int u = current->vertex;
            pq.deleteMin();

            if (visited[u]) continue;

            visited[u] = true;

            for (int i = 0; i < V; i++) {
                const auto& neighbor = adjacencyList[u][i];
                int v = neighbor.first;
                int weight = neighbor.second;

                if (!visited[v] && distance[u] != numeric_limits<int>::max() &&
                    distance[u] + weight < distance[v]) {
                    distance[v] = distance[u] + weight;

                    int index = -1;
                    if (pq.find(v)) {
                        index = v;
                    }

                    if (index != -1) {
                        pq.decreaseKey(index, distance[v]); 
                    }
                    else {
                        pq.build(new Target(v, distance[v]));
                    }
                }
            }
            delete current;
        }


        auto end = chrono::high_resolution_clock::now();

        // Wypisujemy obliczone odległości
        /*cout << "Distances from source vertex " << source << endl;
        for (int i = 0; i < V; ++i) {
            if (distance[i] == numeric_limits<int>::max()) {
                cout << "Vertex " << i << ": Unreachable" << endl;
            }
            else {
                cout << "Vertex " << i << ": " << distance[i] << endl;
            }
        }*/

        delete[] distance;
        delete[] visited;
        return chrono::duration<double, milli>(end - start).count();
    }
};

class Graph {
private:
    int V;
    AdjacencyList adjacencyList;
    AdjacencyMatrix adjacencyMatrix;
public:
    Graph(int vertices) : V(vertices), adjacencyList(vertices), adjacencyMatrix(vertices) {}

    // MATRIX
    void displayMatrix() {
        adjacencyMatrix.display();
    }

    void displayEdgesMatrix(int vertex) {
		adjacencyMatrix.displayEdges(vertex);
	}

    double dijkstraMatrix(int source, ChangeablePriorityQueue& cpq) {
        return adjacencyMatrix.dijkstra(source, cpq);
    }

    // LIST
    void displayList() {
        adjacencyList.display();
    }

    void displayEdgesList(int vertex) {
        adjacencyList.displayEdges(vertex);
    }

    double dijkstraList(int source, ChangeablePriorityQueue& cpq) {
        return adjacencyList.dijkstra(source, cpq);
    }

    // BOTH
    void display() {
        cout << "Matrix\n";
        displayMatrix();
        cout << "List\n";
		displayList();
    }

    void displayEdges(int vertex) {
        cout << "Edges from vertex: " << vertex << endl;
        cout << "Matrix\n";
		displayEdgesMatrix(vertex);
        cout << "List\n";
		displayEdgesList(vertex);
	}

    double* dijkstra(int source, ChangeablePriorityQueue& cpq) {
        double* distances = new double[2];
		distances[0] = adjacencyMatrix.dijkstra(source, cpq);
		distances[1] = adjacencyList.dijkstra(source, cpq);
        return distances;
	}

    void setVertexLabel(int vertex, string label) {
        adjacencyMatrix.setVertexLabel(vertex, label);
        adjacencyList.setVertexLabel(vertex, label);
    }

    void addRandomEdges(int density, int maxWeight) {
        srand(time(0));
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                if (i != j) { // bo d(s,s) jest rowne 0
                    if (rand() % 100 < density) {
                        int weight = rand() % maxWeight + 1;
                        adjacencyList.addEdge(i, j, weight);
                        adjacencyList.addEdge(j, i, weight);
                        adjacencyMatrix.addEdge(i, j, weight);
                        adjacencyMatrix.addEdge(j, i, weight);
                    }
                }
            }
        }
    }

};

int main() {
    int numInstances = 1;
    int maxWeight = 10;

    for (int i = 0; i < numInstances; i++) {
		Graph g(5);
		g.setVertexLabel(0, "A");
		g.setVertexLabel(1, "B");
		g.setVertexLabel(2, "C");
		g.setVertexLabel(3, "D");
		g.setVertexLabel(4, "E");

		g.addRandomEdges(75, maxWeight);
		g.display();

		ChangeablePriorityQueue cpq;
        g.displayEdges(0);
        double* time = g.dijkstra(0, cpq);
        cout << "Dijkstra - Source: 0\n";
        cout << "Matrix: " << *(time) << endl;
		cout << "List: " << *(time+1) << endl;
		cout << endl;
        delete[] time;
	}

    getchar();
    return 0;
}

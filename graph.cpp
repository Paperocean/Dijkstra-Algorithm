#include <vector>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <iomanip>

#include <fstream>


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
    int V;

    AdjacencyMatrix(int numVertices) : V(numVertices) {
        adjacencyMatrix = new int* [V];
        for (int i = 0; i < V; ++i) {
            adjacencyMatrix[i] = new int[V];
            for (int j = 0; j < V; ++j) {
                adjacencyMatrix[i][j] = -1;
            }
        }
    }

    ~AdjacencyMatrix() {
        for (int i = 0; i < V; ++i) {
            delete[] adjacencyMatrix[i];
        }
        delete[] adjacencyMatrix;
    }

    void addEdge(int u, int v, int w) {
        adjacencyMatrix[u][v] = w;
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

            delete current; // Zwolnienie pamięci zaalokowanej dla wierzchołka
        }

        auto end = chrono::high_resolution_clock::now();
        return chrono::duration<double, milli>(end - start).count();

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
    }
};
struct AdjacencyList {
    pair<int, int>** adjacencyList;
    int V;

    AdjacencyList(int numVertices) : V(numVertices) {
        adjacencyList = new pair<int, int>* [V];
        for (int i = 0; i < V; i++) {
            adjacencyList[i] = new pair<int, int>[V];
            for (int j = 0; j < V; j++) {
                adjacencyList[i][j] = make_pair(-1, -1);
            }
        }
    }

    ~AdjacencyList() {
        for (int i = 0; i < V; i++) {
            delete[] adjacencyList[i];
        }
        delete[] adjacencyList;
    }

    void addEdge(int u, int v, int w) {
        adjacencyList[u][v] = make_pair(v, w);
        adjacencyList[v][u] = make_pair(u, w);
    }

    void display() {
        for (int i = 0; i < V; i++) {
            cout << "Vertex " << i << ": ";
            for (int j = 0; j < V; j++) {
                const auto& neighbor = adjacencyList[i][j];
                cout << "(" << neighbor.first << ", " << neighbor.second << ") ";
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

        int* distance = new int[V]; // Tablica wag/odleglosci dla wierzcholkow
        bool* visited = new bool[V]; // Tablica odwiedzonych wierzchołków

        for (int i = 0; i < V; ++i) {
            distance[i] = numeric_limits<int>::max();
            visited[i] = false;
        }
        distance[source] = 0; // Odległość od źródła do źródła wynosi 0

        pq.build(new Target(source, 0)); // Umieszczamy źródło w kolejce

        // Algorytm Dijkstry
        while (!pq.isEmpty()) {
            Target* current = pq.peek(); // Wybieramy wierzchołek o najmniejszej odległości
            int u = current->vertex;
            pq.deleteMin(); // Usuwamy wybrany wierzchołek z kolejki

            if (visited[u]) continue; // Jeśli wierzchołek był już odwiedzony, przechodzimy do następnego

            visited[u] = true; // Oznaczamy wierzchołek jako odwiedzony

            // Aktualizujemy odległości dla sąsiadów wybranego wierzchołka
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
                        pq.decreaseKey(index, distance[v]); // Aktualizujemy odległość dla v
                        // Dlaczego? Bo waga krawedzi była wieksza 
                        // niz dodanie do siebie wagi od wierzcholka zrodlowego
                    }
                    else {
                        pq.build(new Target(v, distance[v])); // Dodajemy v do kolejki
                        // Dodajesz nowe v do kolejki o tym dystansie 
                        // bo nie znalazles go w kolejce
                    }
                }
            }
            delete current; // Zwolnienie pamięci zaalokowanej dla wierzchołka
        }


        auto end = chrono::high_resolution_clock::now();
        return chrono::duration<double, milli>(end - start).count();

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
    }
};

class Graph {
private:
    int V;
    AdjacencyList adjacencyList;
    AdjacencyMatrix adjacencyMatrix;
public:
    Graph(int vertices) : V(vertices), adjacencyList(vertices), adjacencyMatrix(vertices) {}

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

    void displayMatrix() {
        adjacencyMatrix.display();
    }

    void displayList() {
        adjacencyList.display();
    }

    double dijkstraMatrix(int source, ChangeablePriorityQueue& cpq) {
        return adjacencyMatrix.dijkstra(source, cpq);
    }

    double dijkstraList(int source, ChangeablePriorityQueue& cpq) {
        return adjacencyList.dijkstra(source, cpq);
    }

};

int main() {
    int numInstances = 1; // Zmiana liczby instancji na 100
    int maxWeight = 10; // maksymalna waga krawedzi
    int vMax = 100; // maksymalna ilosc wierzcholkow

    vector<int> numVerticesList = { 5, 10, 15, 20, 25 };
    vector<int> densities = { 25, 50, 75, 100 };
    double timeTakenM{ 0 }, timeTakenL{ 0 };

    for (int numVertices : numVerticesList) {
        double timeTakenAvgM{ 0 }, timeTakenAvgL{ 0 };
        for (int density : densities) {
            cout << "Adjacency Graph - Num Vertices: " << numVertices << ", Density: " << density << endl;
            for (int i = 0; i < numInstances; i++) {
                Graph g(numVertices);
                g.addRandomEdges(density, maxWeight);
                g.displayMatrix();
                g.displayList(); 

                ChangeablePriorityQueue cpq;

                //timeTakenM = g.dijkstraMatrix(0, cpq);
                //timeTakenL = g.dijkstraList(0, cpq);

                //// Wyswietlenie tego co zostaje zapisane do csv
                //// Wyswietlenie czasow dla KAZDEGO pomiaru 
                //// dla jednego wierzcholka zrodlowego
                //cout << "Matrix: " << timeTakenM << endl;
                //cout << "List: " << timeTakenL << endl;
                //cout << endl;

                // Wyswietlenie czasow dla KAZDEGO pomiaru 
                // dla wszystkich mozliwych wierzcholkow zrodlowych
                for (int j = 0; j < numVertices; j++) {
                    timeTakenM = g.dijkstraMatrix(j, cpq);
                    timeTakenL = g.dijkstraList(j, cpq);
                    cout << "Dijkstra - Source: " << j << " Vertices: " << numVertices << " Density: " << density << endl;
                    cout << "Matrix: " << timeTakenM << endl;
                    cout << "List: " << timeTakenL << endl;
                    cout << endl;
                }

                // czasy usrednione
                timeTakenAvgM += timeTakenM;
                timeTakenAvgL += timeTakenL;
            }
        }
        // Wyswietlenie usrednionych wartosci pod sprawozdanie
        cout << "Average time taken\nMatrix: " << fixed << setprecision(4) << timeTakenAvgM / numInstances << " ms\nList: " << timeTakenAvgL / numInstances << " ms" << endl;
    }
    getchar();
    return 0;
}
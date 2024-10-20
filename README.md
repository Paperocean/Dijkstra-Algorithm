# Dijkstra's Algorithm

## Overview
Dijkstra's algorithm is a popular algorithm used for finding the shortest paths between nodes in a graph. It is particularly effective for graphs with non-negative edge weights. This repository contains an implementation of Dijkstra's algorithm, which can be used to solve various pathfinding problems in weighted graphs.

## Features
- **Shortest Path Calculation:** Efficiently calculates the shortest path from a source node to all other nodes in the graph.
- **Non-Negative Edge Weights:** Designed to work with graphs that have non-negative weights, ensuring accurate results.
- **Graph Representation:** Supports graph representation using adjacency lists or matrices.

## Time Complexity
- **Time Complexity:** O(V^2) with a simple implementation, where V is the number of vertices. Using a priority queue, the complexity can be improved to O(E + V log V), where E is the number of edges.

## Usage
To use the Dijkstra's algorithm implementation in this repository:
1. Clone the repository:
   ```bash
   git clone https://github.com/Paperocean/Dijkstra-algorithm.git
2. Navigate to the project directory:
   ```bash
   cd Dijkstra-algorithm
3. Run the algorithm with your graph data. Refer to the provided examples for guidance.

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

public class shortestPathAlgo {

  private final int N, start;
  private final double[][] distance;
  private List<Integer> tour = new ArrayList<>();
  private double shortestDistance = Double.POSITIVE_INFINITY;
  private boolean status = false;

  public shortestPathAlgo(double[][] distance) 
  {
    this(0, distance);
  } 

  public shortestPathAlgo(int start, double[][] distance) 
  {
    N = distance.length;
    
//error checking
    if (N <= 2) throw new IllegalStateException("N <= 2! Invalid.");
    if (N != distance[0].length) throw new IllegalStateException("Matrix should be a square matrix");
    if (start < 0 || start >= N) throw new IllegalArgumentException("Invalid Start Node");

    this.start = start;
    this.distance = distance;
  }

  
// Returns the best route for the tour
  public List<Integer> getPath() {
    if (!status) solve();
    return tour;
  }

  // Returns the shortest tour path
  public double getDistance() 
  {
    if (!status) solve();
    return shortestDistance;
  }

  // Algorithm to find the shortest tour
  
  public void solve() {

      
    if (status) return;

    
    final int END_STATE = (1 << N) - 1;
    
    Double[][] memo = new Double[N][1 << N];

    //Adding nodes to a given node
    
    for (int end = 0; end < N; end++) 
    {
      if (end == start) continue;
      memo[end][(1 << start) | (1 << end)] = distance[start][end];
    }

    for (int r = 3; r <= N; r++) 
    {
      for (int subset : combinations(r, N)) {
          
        if (notIn(start, subset)) continue;
        for (int next = 0; next < N; next++) 
        {
          if (next == start || notIn(next, subset)) continue;
          int subsetWithoutNext = subset ^ (1 << next);
          double minDist = Double.POSITIVE_INFINITY;
          for (int end = 0; end < N; end++)
          {
              
            if (end == start || end == next || notIn(end, subset)) continue;
            
            double newDistance = memo[end][subsetWithoutNext] + distance[end][next];
            
            if (newDistance < minDist) 
            {
              minDist = newDistance;
            }
          }
          memo[next][subset] = minDist;
        }
      }
    }

    // Connects the path to strating node back 
    for (int i = 0; i < N; i++) 
    {
      if (i == start) continue;
      
      double tourDistance = memo[i][END_STATE] + distance[i][start];
      
      if (tourDistance < shortestDistance) 
      {
        shortestDistance = tourDistance;
      }
    }

    int lastIndex = start;
    int state = END_STATE;
    tour.add(start);

    // Reconstruct TSP path from memo table.
    for (int i = 1; i < N; i++) 
    {

      int index = -1;
      for (int j = 0; j < N; j++) 
      {
        if (j == start || notIn(j, state)) continue;
        if (index == -1) index = j;
        double prevDist = memo[index][state] + distance[index][lastIndex];
        double newDist  = memo[j][state] + distance[j][lastIndex];
        
        if (newDist < prevDist) 
        {
          index = j;
        }
      }

      tour.add(index);
      state = state ^ (1 << index);
      lastIndex = index;
    }

    tour.add(start);
    Collections.reverse(tour);

    status = true;
  }

  private static boolean notIn(int elem, int subset) {
    return ((1 << elem) & subset) == 0;
  }

 
  public static List<Integer> combinations(int r, int n) {
    List<Integer> subsets = new ArrayList<>();
    combinations(0, 0, r, n, subsets);
    return subsets;
  }

  // To find all the combinations of size r we need to recurse until we have
  // selected r elements (aka r = 0), otherwise if r != 0 then we still need to select
  // an element which is found after the position of our last selected element
  private static void combinations(int set, int at, int r, int n, List<Integer> subsets) {

    // Return early if there are more elements left to select than what is available.
    int elementsLeftToPick = n - at;
    if (elementsLeftToPick < r) return;

    // We selected 'r' elements so we found a valid subset!
    if (r == 0) {
      subsets.add(set);
    } else {
      for (int i = at; i < n; i++) {
        // Try including this element
        set |= 1 << i;

        combinations(set, i + 1, r - 1, n, subsets);

        // Backtrack and try the instance where we did not include this element
        set &= ~(1 << i);
      }
    }
  }

  public static void main(String[] args) {
    

     Scanner myObj = new Scanner(System.in);  // Create a Scanner object
     System.out.println("Enter Distancee between K and L :");
    int distance = myObj.nextInt(); 
    int n = 12;
    double[][] distanceMatrix = new double[n][n];
    for (double[] row : distanceMatrix) java.util.Arrays.fill(row, 10000);
    
    //creating the adjascency matrix
    distanceMatrix[10][11] =distance ;
    distanceMatrix[11][10] =distance ;
    
 //filling the adjecency matrix
    distanceMatrix[0][1] =18 ;
    distanceMatrix[1][0] =18 ;
    distanceMatrix[0][7] =15 ;
    distanceMatrix[7][0] =15 ;
    distanceMatrix[0][8] =20 ;
    distanceMatrix[8][0] =20 ;
    distanceMatrix[1][0] =18 ;
    distanceMatrix[0][1] =18 ;
    distanceMatrix[1][2] = 30;
    distanceMatrix[2][1] = 30;
    distanceMatrix[1][6] = 14;
    distanceMatrix[6][1] = 14;
    distanceMatrix[1][7] = 12;
    distanceMatrix[7][1] = 12;
    distanceMatrix[2][1] = 30;
    distanceMatrix[1][2] = 30;
    distanceMatrix[2][3] = 24;
    distanceMatrix[3][2] = 24;
    distanceMatrix[2][4] = 33;
    distanceMatrix[4][4] = 33;
    distanceMatrix[2][5] = 23;
    distanceMatrix[2][11] = 35;
    distanceMatrix[3][2] = 24;
    distanceMatrix[3][4] = 27;
    distanceMatrix[4][2] = 33;
    distanceMatrix[4][3] = 27;
    distanceMatrix[4][11] = 9;
    distanceMatrix[5][2] = 23;
    distanceMatrix[5][6] = 16;
    distanceMatrix[5][10] = 10;
    distanceMatrix[5][11] = 17;
    distanceMatrix[6][1] = 14;
    distanceMatrix[6][5] = 16;
    distanceMatrix[6][9] = 10;
    distanceMatrix[7][0] = 15;
    distanceMatrix[7][1] = 12;
    distanceMatrix[7][8] = 13;
    distanceMatrix[7][9] = 11;
    distanceMatrix[8][0] = 20;
    distanceMatrix[8][7] = 13;
    distanceMatrix[8][9] = 10;
    distanceMatrix[9][6] = 10;
    distanceMatrix[9][7] = 11;
    distanceMatrix[9][8] = 10;
    distanceMatrix[9][10] = 23;
    distanceMatrix[10][5] = 22;
    distanceMatrix[10][9] = 23;
    distanceMatrix[11][2] = 35;
    distanceMatrix[11][4] = 9;
    distanceMatrix[11][5] = 17;
   

    
    int startNode;
    char userInput;
    Scanner s=new Scanner(System.in);
    System.out.println("Enter starting city : ");
    char city=s.next().charAt(0);
    startNode=city;
    startNode=startNode-65;
    shortestPathAlgo solver = new shortestPathAlgo(startNode, distanceMatrix);

//priting the output
    ArrayList<Integer> temp=new ArrayList<Integer> (solver.getPath());
    for(int i=0;i<temp.size();i++)
    {
        int t=temp.get(i);
        t=t+65;
       char c=(char)t;
       if(t!=startNode+65 || i==0 )
       {
        System.out.print(c+"->");
       }
       else
           System.out.print(c);
          
    
    }

    // Print: 42.0
    System.out.println("\nShortest Path Distance: " + solver.getDistance());
  }
}
import java.util.Random;

/**
 * Created with IntelliJ IDEA.
 * User: Gavin Gresham
 * Date: 9/9/13
 * Time: 12:30 AM
 */
public class RandomSequence implements ArrayGenerator
{
    /**
     *
     * @param N The 2^N power of the array
     * @return A sequential array of length 2^N
     */
    @Override
    public int[] getArray(int N)
    {
        int elements = (int)Math.pow(2,N);

        int[] retn = new int[elements];

        //Fisher-Yates Shuffle
        for(int i = 0; i < elements - 1; i++)
            retn[i] = i + 1;

        //Set the end of the array
        retn[elements - 1] = -1;

        Random rand = new Random();
        //To initialize an array a of n elements to a randomly shuffled copy of source, both 0-based:
        for(int i = 1; i < elements; i++)
        {
            int num = rand.nextInt(i + 1);
            if(num != i)
            {
                int temp = retn[i];
                retn[i] = retn[num];
                retn[num] = temp;
            }
        }

        //Attach any disconnected elements
        int count = 0;
        while(count != elements)
        {
            char[] visited = new char[elements];
            int lastidx = 0;
            int idx = 0;
            count = 0;
            while(idx != -1 && visited[idx] == 0)
            {
                count++;
                visited[idx] = 1;
                lastidx = idx;
                idx = retn[idx];
            }

            //Unroll cycles
            if(count != elements)
            {
                //Find first unvisited
                int visitedIdx = 0;
                while(visited[++visitedIdx] != 0);

                int nextIdx = retn[visitedIdx];
                retn[lastidx] = nextIdx;
                retn[visitedIdx] = -1;

            }
            else
                retn[lastidx] = -1;
        }

        return retn;
    }

    @Override
    public String toString()
    {
        return "Random";
    }
}

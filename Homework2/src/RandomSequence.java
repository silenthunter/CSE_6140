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

        return retn;
    }

    @Override
    public String toString()
    {
        return "Random";
    }
}

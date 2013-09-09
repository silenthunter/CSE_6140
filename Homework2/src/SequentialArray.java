/**
 * Created with IntelliJ IDEA.
 * User: Gavin Gresham
 * Date: 9/9/13
 * Time: 12:35 AM
 */
public class SequentialArray implements ArrayGenerator
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

        for(int i = 0; i < elements; i++)
        {
            retn[i] = i + 1;
        }

        retn[retn.length -1] = -1;

        return retn;
    }

    @Override
    public String toString()
    {
        return "Sequential";
    }
}

/**
 * Created with IntelliJ IDEA.
 * User: Gavin Gresham
 * Date: 9/9/13
 * Time: 12:50 AM
 */
public class NStride implements ArrayGenerator
{
    private int stride;

    public NStride(int stride)
    {
        this.stride = stride;
    }

    @Override
    public int[] getArray(int N)
    {
        int elements = (int)Math.pow(2, N);

        int[] retn = new int[elements];

        for(int i = 0 ; i < elements; i++)
        {
            int val = (i + stride - 1) % (elements - 1) + 1;
            retn[i] = val;
        }

        retn[elements - 1] = -1;

        return retn;
    }

    @Override
    public String toString()
    {
        return "Stride" + stride;
    }
}

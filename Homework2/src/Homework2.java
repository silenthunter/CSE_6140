import java.util.Calendar;
import java.util.Date;
import java.util.ArrayList;
import java.util.Random;

/**
 * Created with IntelliJ IDEA.
 * User: Gavin Gresham
 * Date: 9/8/13
 * Time: 2:06 PM
 */
public class Homework2
{

    public long traverseArray(int arr[])
    {
        long total = 0;
        int idx = 0;
        while((idx = arr[idx]) != -1)
        {
            total += idx;
        }

        return total;
    }

    public long sumOfN(long N)
    {
        long retn = N * (N + 1) / 2;

        return retn;
    }

    public void runSequence(ArrayGenerator arrayClass) throws Exception
    {
        final int average = 5;

        for(int i = 5; i <= 25; i++)
        {
            long totalDiff = 0;
            for(int j = 0; j < average; j++)
            {
                int[] arr = arrayClass.getArray(i);

                long startTime = System.currentTimeMillis();
                long sum = traverseArray(arr);
                long endTime = System.currentTimeMillis();

                //Sum of 1 - (N - 1), because the last element of the array will be -1
                long sumOfN = sumOfN((long) Math.pow(2, i) - 1);

                if(sum != sumOfN)
                    throw new Exception();

                long diff = endTime - startTime;
                totalDiff += diff;
                System.gc();
            }
            totalDiff /= average;
            long onlySeconds = totalDiff % 60000;
            long milli = onlySeconds % 1000;
            long seconds = onlySeconds / 1000;

            System.out.print(arrayClass.toString() + ",");
            System.out.print(i);
            System.out.println("," + onlySeconds);
            //System.out.println(" Milli: " + milli);
        }
    }

    public void timeLinkList()
    {
        final int average = 5;

        for(int i = 5; i <= 25; i++)
        {
            long totalDiff = 0;
            for(int j = 0; j < average; j++)
            {
                int elements = (int)Math.pow(2, i);
                Node root = new Node();
                Node cur = root;
                for(int k = 0; k < elements; k++)
                {
                    cur.next = new Node();
                    cur = cur.next;
                }

                cur = root;
                long startTime = System.currentTimeMillis();
                while((cur = cur.next) != null);
                long endTime = System.currentTimeMillis();


                long diff = endTime - startTime;
                totalDiff += diff;
                System.gc();
            }
            totalDiff /= average;
            long onlySeconds = totalDiff % 60000;
            long milli = onlySeconds % 1000;
            long seconds = onlySeconds / 1000;

            System.out.print("LinkedList,");
            System.out.print(i);
            System.out.println("," + onlySeconds);
            //System.out.println("," + milli);
        }
    }

    class Node
    {
        public Node next = null;
    }

    public static void main(String args[])
    {
        Homework2 hw = new Homework2();

        try
        {
            //hw.timeLinkList();
            //hw.runSequence(new SequentialArray());
            //hw.runSequence(new RandomSequence());
            hw.runSequence(new NStride(2));
            hw.runSequence(new NStride(4));
        }catch(Exception e)
        {
            e.printStackTrace();
        }
    }
}

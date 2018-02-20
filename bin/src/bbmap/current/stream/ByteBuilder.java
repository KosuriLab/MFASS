package stream;

import java.util.Arrays;

import align2.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 8, 2013
 *
 */
public final class ByteBuilder {
	
	public static void main(String[] args){
		StringBuilder sb=new StringBuilder();
	}
	
	public ByteBuilder(){
		array=new byte[32];
	}
	
	public ByteBuilder(int initial){
		assert(initial>=1);
		array=new byte[initial];
	}
	
	public ByteBuilder(Object o){
		String s=o.toString();
		array=new byte[s.length()+1];
		append(s);
	}


	public ByteBuilder append(float x, int places){return append(String.format("%."+places+"f", x));}
	public ByteBuilder append(double x, int places){return append(String.format("%."+places+"f", x));}
	
	public ByteBuilder append(float x){return append(Float.toString(x));}
	public ByteBuilder append(double x){return append(Double.toString(x));}
	public ByteBuilder append(boolean x){return append(x ? tbool : fbool);}
	
	
	public ByteBuilder append(char x){
		if(length>=array.length){expand();}
		array[length]=(byte)x;
		length++;
		return this;
	}
	public ByteBuilder append(byte x){
		if(length>=array.length){expand();}
		array[length]=x;
		length++;
		return this;
	}
	
	public ByteBuilder append(int x){
		expand(11);
		if(x<0){
			if(x==Integer.MIN_VALUE){
				return append(Integer.toString(Integer.MIN_VALUE));
			}else{
				array[length]='-';
				length++;
				x=-x;
			}
		}else if(x==0){
			array[length]='0';
			length++;
			return this;
		}

//		final int len=lengthOf(x);
//		int pos=length+len-1;
//		while(x>9){
//			int y=x%100;
//			x=x/100;
//			array[pos]=ones100[y];
//			pos--;
//			array[pos]=tens100[y];
//			pos--;
//		}
//		while(x>0){
//			int y=x%10;
//			x=x/10;
//			array[pos]=numbers[y];
//			pos--;
//		}
//		length+=len;
		
//		final int initial=length;
//		while(x>9){
//			int y=x%100;
//			x=x/100;
//			array[length]=tens100[y];
//			length--;
//			array[length]=ones100[y];
//			length--;
//		}
//		while(x>0){
//			int y=x%10;
//			x=x/10;
//			array[length]=numbers[y];
//			length++;
//		}
//		
//		for(int i=initial, j=length-1; i<j; i++, j--){
//			byte temp=array[i];
//			array[i]=array[j];
//			array[j]=temp;
//		}
		

		
		int pos=0;
		while(x>9){
			int y=x%100;
			x=x/100;
			numbuffer[pos]=ones100[y];
			pos++; 
			numbuffer[pos]=tens100[y];
			pos++;
		}
		while(x>0){
			int y=x%10;
			x=x/10;
			numbuffer[pos]=ones100[y];
			pos++;
		}
		
		while(pos>0){
			pos--;
			array[length]=numbuffer[pos];
			length++;
		}
		
		return this;
	}
	
	public ByteBuilder append(long x){
		if(x>Integer.MIN_VALUE && x<=Integer.MAX_VALUE){return append((int)x);}
		expand(20);
		if(x<0){
			if(x==Integer.MIN_VALUE){
				return append((long)x);
			}else{
				array[length]='-';
				length++;
				x=-x;
			}
		}else if(x==0){
			array[length]='0';
			length++;
			return this;
		}

//		final int len=lengthOf(x);
//		int pos=length+len-1;
//		while(x>9){
//			int y=(int)(x%100);
//			x=x/100;
//			array[pos]=ones100[y];
//			pos--;
//			array[pos]=tens100[y];
//			pos--;
//		}
//		while(x>0){
//			int y=(int)(x%10);
//			x=x/10;
//			array[pos]=numbers[y];
//			pos--;
//		}
//		length+=len;
		
		int pos=0;
		while(x>9){
			int y=(int)(x%100);
			x=x/100;
			numbuffer[pos]=ones100[y];
			pos++;
			numbuffer[pos]=tens100[y];
			pos++;
		}
		while(x>0){
			int y=(int)(x%10);
			x=x/10;
			numbuffer[pos]=ones100[y];
			pos++;
		}
		
		while(pos>0){
			pos--;
			array[length]=numbuffer[pos];
			length++;
		}
		
		return this;
	}
	
	public ByteBuilder append(String x){
		if(x==null){return append(nullBytes);}
		expand(x.length());
		for(int i=0; i<x.length(); i++){
			array[length]=(byte)x.charAt(i);
			length++;
		}
		return this;
	}
	
	public ByteBuilder append(CharSequence x){
		if(x==null){return append(nullBytes);}
		expand(x.length());
		for(int i=0; i<x.length(); i++){
			array[length]=(byte)x.charAt(i);
			length++;
		}
		return this;
	}
	
	public ByteBuilder append(Object x){
		if(x==null){return append(nullBytes);}
		return append(x.toString());
	}
	
	public ByteBuilder append(byte[] x){
		if(x==null){x=nullBytes;}
		expand(x.length);
		for(int i=0; i<x.length; i++){
			array[length]=x[i];
			length++;
		}
		return this;
	}
	
	public ByteBuilder append(byte[] x, int start, int len){
//		if(x==null){x=nullBytes;}
		expand(len);
		for(int i=start, lim=start+len; i<lim; i++){
			array[length]=x[i];
			length++;
		}
		return this;
	}
	
	public ByteBuilder append(char[] x){
		if(x==null){return append(nullBytes);}
		expand(x.length);
		for(int i=0; i<x.length; i++){
			array[length]=(byte)x[i];
			length++;
		}
		return this;
	}
	
	@Override
	public final String toString(){
		return new String(array, 0, length);
	}
	
	private final boolean isRoom(int x){
		return array.length-length>=x;
	}
	
	private final void expand(){
		int x=Tools.min(Integer.MAX_VALUE, array.length*2);
		array=Arrays.copyOf(array, x);
	}
	
	private final void expand(int extra){
		long x=array.length;
		while(x-length<extra){x<<=1;}
		if(x!=array.length){
			array=Arrays.copyOf(array, (int)Tools.min(Integer.MAX_VALUE, x));
		}
	}
	
	public final void ensureExtra(int extra){
		if(array.length-length<extra){expand(extra);}
	}

	public int length(){return length;}
	public void setLength(int x){
		assert(x>=0 && x<=array.length);
		length=x;
	}

	public static final byte[] numbers=new byte[] {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
	public static final byte[] nullBytes="null".getBytes();
	public static final byte[] fbool="false".getBytes();
	public static final byte[] tbool="true".getBytes();
	public byte[] array;
	public int length=0;
	private final byte[] numbuffer=new byte[19];

	public static final byte[] ones100, tens100;
	
	static{
		ones100=new byte[100];
		tens100=new byte[100];
		for(int i=0; i<100; i++){
			ones100[i]=(byte)('0'+i%10);
			tens100[i]=(byte)('0'+i/10);
		}
	}
	
}

package at.cibiv.argos.io;

import java.io.IOException;
import java.io.InputStream;

/**
 * A stream of items that supports the pop() operation.
 * @author niko.popitsch@univie.ac.at
 * 
 * @param <T>
 */
public interface PopableStream<T> {

	/**
	 * Pop data from the stream.
	 * 
	 * @return
	 */
	public T pop() throws IOException;

	/**
	 * Close the stream.
	 * 
	 * @throws IOException
	 */
	public void close() throws IOException;
	
	/**
	 * @return the number of available bytes.
	 * @throws IOException
	 */
	public int available() throws IOException;
	
	/**
	 * @return the stream data type
	 */
	public Class<T> getDatatype();
	
	/**
	 * @return the number of decoded entities if available
	 */
	public Long getDecodedEntities();
	
	/**
	 * @return the underlying input stream if any.
	 */
	public InputStream getInputStream();

	/**
	 * 
	 * @return the number of bytes read from this stream
	 */
	public Long getReadBytes();
}

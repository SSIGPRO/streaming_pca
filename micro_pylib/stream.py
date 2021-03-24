import serial
import struct
import numpy as np

class Stream:
    def __init__(self, com_port, baud_rate=115200, timeout=1, double=False):
        print("Opening serial port to micro...")
        self.serial_port = serial.Serial(
            com_port,
            baud_rate,
            bytesize=serial.EIGHTBITS,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE,
            timeout=timeout,
            )
        self.double = double
    
    def __del__(self):
        self.close()
        
    def close(self):
        print("Closing serial port to micro...")
        self.serial_port.close()
    
    def write_matrix(self, M, send_shape=False):
        if send_shape:
            # write the shape of the matrix
            self.write_value(M.shape[0])
            self.write_value(M.shape[1])
        # Pack
        M_list = list(M.flatten())
        if self.double:
            bytes = struct.pack("%sd"%(M.shape[0]*M.shape[1]), *M_list)
        else:
            bytes = struct.pack("%sf"%(M.shape[0]*M.shape[1]), *M_list)
        # write
        self.serial_port.write(bytes)
        
    def write_vector(self, M, send_shape=False):
        if send_shape:
            # write the shape of the matrix
            self.write_value(M.shape[0])
        # Pack
        M_list = list(M.flatten())
        if self.double:
            bytes = struct.pack("%sd"%(M.shape[0]), *M_list)
        else:
            bytes = struct.pack("%sf"%(M.shape[0]), *M_list)
        # write
        self.serial_port.write(bytes)
        
    def write_value(self, val, type="int"):
        # Pack
        bytes = None
        if type == "int":
            bytes = struct.pack("I", val)
        else:
            if self.double:
                bytes = struct.pack("d", val)
            else:
                bytes = struct.pack("f", val)
        # write
        self.serial_port.write(bytes)

    
    def read_matrix(self, rows=None, cols=None):
        if cols is None or rows is None:
            # Receive the size of the matrix
            n = self.read_value()
            m = self.read_value()
        else:
            n = rows
            m = cols
        length = n*m
        # Receive
        total_bytes = length*(4+4*self.double)
        received_bytes = 0
        bytes = b""
        while received_bytes < total_bytes:
            current_bytes = total_bytes-received_bytes
            if current_bytes > 0xFFFF:
                current_bytes = 0xFFFF
            bytes += self.serial_port.read(current_bytes)
            received_bytes += 0xFFFF
            
        # Unpack
        if self.double:
            M_list = struct.unpack("%sd"%(length), bytes)
        else:
            M_list = struct.unpack("%sf"%(length), bytes)
        M = np.array(M_list)
        M = np.reshape(M, (n,m))
        return M
        
    def read_vector(self, length=None):
        if length is None:
            # Receive the size of the vector
            n = self.read_value()
        else:
            n = length
        # Receive
        bytes = self.serial_port.read(n*(4+4*self.double))
        # Unpack
        if self.double:
            v_list = struct.unpack("%sd"%n, bytes)
        else:
            v_list = struct.unpack("%sf"%n, bytes)
        v = np.array(v_list)
        return v
        
    def read_value(self, type="int"):
        # Receive
        if type!="int":
            bytes = self.serial_port.read(4+4*self.double)
        else:
            bytes = self.serial_port.read(4)
        if bytes != b'':
            # Unpack
            val = None
            if type == "int":
                val = struct.unpack("I", bytes)[0]
            else:
                if self.double:
                    val = struct.unpack("d", bytes)[0]
                else:
                    val = struct.unpack("f", bytes)[0]
            return val
        else:
            return None
            
    def write_msg(self, bytes):
        return self.serial_port.write(bytes)
    
    def read_msg(self, length):
        return self.serial_port.read(length)
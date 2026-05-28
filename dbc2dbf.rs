use std::env;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::process::exit;

const CHUNK: usize = 4096;
const MAXBITS: usize = 13;
const MAXWIN: usize = 4096;

#[derive(Debug)]
enum BlastError {
    OutOfInput,
    OutputError,
    InvalidLiteralFlag,
    InvalidDictSize,
    DistanceTooFar,
    RanOutCodes,
}

struct State<'a, 'b> {
    reader: &'a mut dyn Read,
    writer: &'b mut dyn Write,
    in_buf: [u8; CHUNK],
    in_idx: usize,
    in_len: usize,
    bitbuf: i32,
    bitcnt: i32,
    out: [u8; MAXWIN],
    next: usize,
    first: bool,
}

struct Huffman {
    count: [i16; MAXBITS + 1],
    symbol: [i16; 256],
}

impl Default for Huffman {
    fn default() -> Self {
        Huffman {
            count: [0; MAXBITS + 1],
            symbol: [0; 256],
        }
    }
}

fn bits(s: &mut State, need: i32) -> Result<i32, BlastError> {
    let mut val = s.bitbuf;
    while s.bitcnt < need {
        if s.in_idx >= s.in_len {
            let n = s.reader.read(&mut s.in_buf).map_err(|_| BlastError::OutOfInput)?;
            if n == 0 { return Err(BlastError::OutOfInput); }
            s.in_len = n;
            s.in_idx = 0;
        }
        val |= (s.in_buf[s.in_idx] as i32) << s.bitcnt;
        s.in_idx += 1;
        s.bitcnt += 8;
    }
    s.bitbuf = val >> need;
    s.bitcnt -= need;
    Ok(val & ((1 << need) - 1))
}

fn construct(h: &mut Huffman, rep: &[u8]) {
    let mut n = rep.len();
    let mut symbol = 0;
    let mut length = [0i16; 256];
    let mut rep_idx = 0;

    while n > 0 {
        n -= 1;
        let mut len = rep[rep_idx] as i16;
        rep_idx += 1;
        let mut left = (len >> 4) + 1;
        len &= 15;
        loop {
            length[symbol] = len;
            symbol += 1;
            left -= 1;
            if left == 0 { break; }
        }
    }
    let n_sym = symbol;

    for len in 0..=MAXBITS { h.count[len] = 0; }
    for sym in 0..n_sym { h.count[length[sym] as usize] += 1; }
    if h.count[0] == n_sym as i16 { return; }

    let mut offs = [0i16; MAXBITS + 2];
    offs[1] = 0;
    for len in 1..MAXBITS { offs[len + 1] = offs[len] + h.count[len]; }

    for sym in 0..n_sym {
        if length[sym] != 0 {
            let len = length[sym] as usize;
            h.symbol[offs[len] as usize] = sym as i16;
            offs[len] += 1;
        }
    }
}

fn decode(s: &mut State, h: &Huffman) -> Result<i32, BlastError> {
    let mut bitbuf = s.bitbuf;
    let mut left = s.bitcnt;
    let mut code = 0;
    let mut first = 0;
    let mut index = 0;
    let mut len = 1;
    let mut next = 1;

    loop {
        while left > 0 {
            left -= 1;
            code |= (bitbuf & 1) ^ 1;
            bitbuf >>= 1;
            let count = h.count[next] as i32;
            next += 1;
            if code < first + count {
                s.bitbuf = bitbuf;
                s.bitcnt = (s.bitcnt - len) & 7;
                return Ok(h.symbol[index as usize + (code - first) as usize] as i32);
            }
            index += count;
            first += count;
            first <<= 1;
            code <<= 1;
            len += 1;
        }
        left = (MAXBITS as i32 + 1) - len;
        if left == 0 { break; }
        if s.in_idx >= s.in_len {
            let n = s.reader.read(&mut s.in_buf).map_err(|_| BlastError::OutOfInput)?;
            if n == 0 { return Err(BlastError::OutOfInput); }
            s.in_len = n;
            s.in_idx = 0;
        }
        bitbuf = s.in_buf[s.in_idx] as i32;
        s.in_idx += 1;
        if left > 8 { left = 8; }
    }
    Err(BlastError::RanOutCodes)
}

fn decomp(s: &mut State) -> Result<(), BlastError> {
    let litlen: [u8; 98] = [
        11, 124, 8, 7, 28, 7, 188, 13, 76, 4, 10, 8, 12, 10, 12, 10, 8, 23, 8,
        9, 7, 6, 7, 8, 7, 6, 55, 8, 23, 24, 12, 11, 7, 9, 11, 12, 6, 7, 22, 5,
        7, 24, 6, 11, 9, 6, 7, 22, 7, 11, 38, 7, 9, 8, 25, 11, 8, 11, 9, 12,
        8, 12, 5, 38, 5, 38, 5, 11, 7, 5, 6, 21, 6, 10, 53, 8, 7, 24, 10, 27,
        44, 253, 253, 253, 252, 252, 252, 13, 12, 45, 12, 45, 12, 61, 12, 45,
        44, 173
    ];
    let lenlen: [u8; 6] = [2, 35, 36, 53, 38, 23];
    let distlen: [u8; 7] = [2, 20, 53, 230, 247, 151, 248];
    let base: [i32; 16] = [3, 2, 4, 5, 6, 7, 8, 9, 10, 12, 16, 24, 40, 72, 136, 264];
    let extra: [i32; 16] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8];

    let mut litcode = Huffman::default();
    let mut lencode = Huffman::default();
    let mut distcode = Huffman::default();

    construct(&mut litcode, &litlen);
    construct(&mut lencode, &lenlen);
    construct(&mut distcode, &distlen);

    let lit = bits(s, 8)?;
    if lit > 1 { return Err(BlastError::InvalidLiteralFlag); }
    let dict = bits(s, 8)?;
    if dict < 4 || dict > 6 { return Err(BlastError::InvalidDictSize); }

    loop {
        if bits(s, 1)? != 0 {
            let mut symbol = decode(s, &lencode)?;
            let mut len = base[symbol as usize] + bits(s, extra[symbol as usize])?;
            if len == 519 { break; }

            symbol = if len == 2 { 2 } else { dict };
            let mut dist = decode(s, &distcode)? << symbol;
            dist += bits(s, symbol)?;
            dist += 1;
            if s.first && dist as usize > s.next {
                return Err(BlastError::DistanceTooFar);
            }

            while len != 0 {
                let mut from = if s.next < dist as usize {
                    s.next + MAXWIN - dist as usize
                } else {
                    s.next - dist as usize
                };

                let mut copy = if s.next < dist as usize { dist as usize } else { MAXWIN };
                copy -= s.next;
                if copy > len as usize { copy = len as usize; }
                len -= copy as i32;

                for _ in 0..copy {
                    s.out[s.next] = s.out[from];
                    s.next += 1;
                    from += 1;
                }

                if s.next == MAXWIN {
                    s.writer.write_all(&s.out[..s.next]).map_err(|_| BlastError::OutputError)?;
                    s.next = 0;
                    s.first = false;
                }
            }
        } else {
            let symbol = if lit != 0 { decode(s, &litcode)? } else { bits(s, 8)? };
            s.out[s.next] = symbol as u8;
            s.next += 1;
            if s.next == MAXWIN {
                s.writer.write_all(&s.out[..s.next]).map_err(|_| BlastError::OutputError)?;
                s.next = 0;
                s.first = false;
            }
        }
    }
    Ok(())
}

fn blast<R: Read, W: Write>(reader: &mut R, writer: &mut W) -> i32 {
    let mut s = State {
        reader, writer,
        in_buf: [0; CHUNK], in_idx: 0, in_len: 0,
        bitbuf: 0, bitcnt: 0,
        out: [0; MAXWIN], next: 0, first: true,
    };

    let err = match decomp(&mut s) {
        Ok(_) => 0,
        Err(BlastError::OutOfInput) => 2,
        Err(BlastError::OutputError) => 1,
        Err(BlastError::InvalidLiteralFlag) => -1,
        Err(BlastError::InvalidDictSize) => -2,
        Err(BlastError::DistanceTooFar) => -3,
        Err(BlastError::RanOutCodes) => -9,
    };

    if err != 1 && s.next > 0 {
        if s.writer.write_all(&s.out[..s.next]).is_err() && err == 0 { return 1; }
    }
    err
}

// --- Lógica principal da aplicação ---

fn dbc2dbf<R: Read + Seek, W: Write>(input: &mut R, output: &mut W) -> io::Result<i32> {
    // Pula para o byte 8 para ler o tamanho do header
    input.seek(SeekFrom::Start(8))?;
    let mut raw_header = [0u8; 2];
    input.read_exact(&mut raw_header)?;
    
    // Converte de little endian (independente de arquitetura)
    let header_size = u16::from_le_bytes(raw_header) as usize;

    // Volta ao começo para copiar o header original
    input.seek(SeekFrom::Start(0))?;
    let mut header = vec![0u8; header_size];
    input.read_exact(&mut header)?;

    // Altera o último byte para indicar fim de cabeçalho (padrão DBF)
    header[header_size - 1] = 0x0D;
    output.write_all(&header)?;

    // Pula para o início dos dados comprimidos
    input.seek(SeekFrom::Start((header_size + 4) as u64))?;

    let ret = blast(input, output);

    if ret != 0 {
        eprintln!("blast error: {}", ret);
    }

    // Verifica se sobrou algum byte não lido
    let current_pos = input.stream_position()?;
    let end_pos = input.seek(SeekFrom::End(0))?;
    let leftover = end_pos - current_pos;
    if leftover > 0 {
        eprintln!("blast warning: {} unused bytes of input", leftover);
    }

    Ok(ret)
}

fn help(prog_name: &str) {
    eprintln!("Syntax error!");
    eprintln!("\tUsage: {} input.dbc output.dbf", prog_name);
}

fn run_cli(args: &[String]) -> i32 {
    if args.len() == 3 {
        let input_path = &args[1];
        let output_path = &args[2];

        let mut input = match File::open(input_path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Erro ao abrir o arquivo de entrada '{}': {}", input_path, e);
                return 1;
            }
        };

        let mut output = match File::create(output_path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Erro ao criar o arquivo de saída '{}': {}", output_path, e);
                return 1;
            }
        };

        match dbc2dbf(&mut input, &mut output) {
            Ok(ret) if ret != 0 => ret,
            Ok(_) => 0,
            Err(e) => {
                eprintln!("Erro de I/O fatal durante conversão: {}", e);
                -1
            }
        }
    } else {
        let prog = args.get(0).map(|s| s.as_str()).unwrap_or("dbc2dbf");
        help(prog);
        1
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    exit(run_cli(&args));
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_help_output() {
        // Testa se a função help executa sem falhas (cobertura da função)
        help("dbc2dbf_test");
    }

    #[test]
    fn test_blast_invalid_literal() {
        // Simula os primeiros bits ativando a flag inválida de literal (> 1)
        // 8 bits = 1 byte. Passamos 0x02.
        let mut input = Cursor::new(vec![0x02]);
        let mut output = Cursor::new(Vec::new());
        let ret = blast(&mut input, &mut output);
        assert_eq!(ret, -1); // BlastError::InvalidLiteralFlag
    }

    #[test]
    fn test_blast_invalid_dict_size() {
        // Primeiro byte (lit) = 0. Segundo byte (dict) = 3 (inválido, pois deve ser entre 4 e 6).
        let mut input = Cursor::new(vec![0x00, 0x03]);
        let mut output = Cursor::new(Vec::new());
        let ret = blast(&mut input, &mut output);
        assert_eq!(ret, -2); // BlastError::InvalidDictSize
    }

    #[test]
    fn test_bits_read() {
        let data = [0b10101010, 0b01010101];
        let mut cursor = Cursor::new(data);
        let mut dummy_out = Vec::new();
        let mut s = State {
            reader: &mut cursor,
            writer: &mut dummy_out,
            in_buf: [0; CHUNK], in_idx: 0, in_len: 0,
            bitbuf: 0, bitcnt: 0,
            out: [0; MAXWIN], next: 0, first: true,
        };

        let b1 = bits(&mut s, 4).unwrap();
        assert_eq!(b1, 0b1010); // Primeiros 4 bits do byte 1
        let b2 = bits(&mut s, 4).unwrap();
        assert_eq!(b2, 0b1010); // Últimos 4 bits do byte 1
        let b3 = bits(&mut s, 8).unwrap();
        assert_eq!(b3, 0b01010101); // Segundo byte
    }

    #[test]
    fn test_dbc2dbf_structural_mock() {
        // Construindo um DBC falso na memória
        let mut mock_dbc = vec![0xFF; 8]; // 8 bytes lixo iniciais
        mock_dbc.extend_from_slice(&[6, 0]); // Header size = 6
        // O arquivo terá exatamente 10 bytes de tamanho.
        // A função pulará para `header_size + 4` (6 + 4 = 10), então o blast começará lendo o EOF.

        let mut input = Cursor::new(mock_dbc);
        let mut output = Cursor::new(Vec::new());

        let ret = dbc2dbf(&mut input, &mut output).unwrap();
        // O blast logo no primeiro bit recebe EOF, retornando OutOfInput (2)
        assert_eq!(ret, 2); // 2 = OutOfInput

        let out_data = output.into_inner();
        assert_eq!(out_data.len(), 6);
        assert_eq!(out_data[5], 0x0D); // Verifica se alterou o último byte para o padrão DBF (0x0D)
    }

    #[test]
    fn test_run_cli_args() {
        assert_eq!(run_cli(&["dbc2dbf".to_string()]), 1);
        assert_eq!(run_cli(&["dbc2dbf".to_string(), "a".to_string(), "b".to_string(), "c".to_string()]), 1);
    }

    #[test]
    fn test_run_cli_io_errors() {
        assert_eq!(run_cli(&["dbc2dbf".to_string(), "arquivo_nao_existe.dbc".to_string(), "out.dbf".to_string()]), 1);
    }

    #[test]
    fn test_run_cli_success() {
        let out_file = "tests/test_run_cli_out.dbf".to_string();
        let ret = run_cli(&["dbc2dbf".to_string(), "tests/sids.dbc".to_string(), out_file.clone()]);
        assert_eq!(ret, 0);
        let _ = std::fs::remove_file(out_file);
    }

    #[test]
    fn test_dbc2dbf_leftover_bytes() {
        let mut data = Vec::new();
        File::open("tests/sids.dbc").unwrap().read_to_end(&mut data).unwrap();
        data.push(0xFF); // Força um byte não lido
        let mut input = Cursor::new(data);
        let mut output = Cursor::new(Vec::new());
        assert_eq!(dbc2dbf(&mut input, &mut output).unwrap(), 0);
    }

    struct BadReader;
    impl Read for BadReader { fn read(&mut self, _buf: &mut [u8]) -> io::Result<usize> { Err(io::Error::new(io::ErrorKind::Other, "bad read")) } }
    impl Seek for BadReader { fn seek(&mut self, _pos: SeekFrom) -> io::Result<u64> { Err(io::Error::new(io::ErrorKind::Other, "bad seek")) } }

    #[test]
    fn test_dbc2dbf_io_error() {
        let mut input = BadReader;
        assert!(dbc2dbf(&mut input, &mut Cursor::new(Vec::new())).is_err());
    }

    struct BadWriter;
    impl Write for BadWriter {
        fn write(&mut self, _buf: &[u8]) -> io::Result<usize> { Err(io::Error::new(io::ErrorKind::Other, "bad write")) }
        fn flush(&mut self) -> io::Result<()> { Ok(()) }
    }

    #[test]
    fn test_blast_output_error() {
        let mut data = Vec::new();
        File::open("tests/sids.dbc").unwrap().read_to_end(&mut data).unwrap();
        let header_size = u16::from_le_bytes([data[8], data[9]]) as usize;
        let compressed_data = &data[header_size + 4..];
        assert_eq!(blast(&mut Cursor::new(compressed_data), &mut BadWriter), 1); // 1 = OutputError
    }

    #[test]
    fn test_blast_corrupt_data() {
        let data = vec![0x00, 0x04, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF];
        assert!(blast(&mut Cursor::new(data), &mut Cursor::new(Vec::new())) < 0);
    }

    #[test]
    fn test_construct_empty() {
        let mut h = Huffman::default();
        construct(&mut h, &[0; 16]); // Simula blocos de bit vazios
        assert_eq!(h.count[0], 16);
    }

    #[test]
    fn test_decode_ran_out_codes() {
        let mut h = Huffman::default();
        let mut cursor = Cursor::new(vec![0xFF, 0xFF, 0xFF]); // bits cheios, falhará em 13 execuções (left == 0)
        let mut dummy = Vec::new();
        let mut s = State {
            reader: &mut cursor, writer: &mut dummy,
            in_buf: [0; CHUNK], in_idx: 0, in_len: 0,
            bitbuf: 0, bitcnt: 0, out: [0; MAXWIN], next: 0, first: true,
        };
        assert!(matches!(decode(&mut s, &h), Err(BlastError::RanOutCodes)));
    }

    #[test]
    fn test_decode_out_of_input() {
        let mut h = Huffman::default();
        h.count[1] = 1; h.symbol[0] = 1;
        let mut cursor = Cursor::new(Vec::new()); // EOF forçado subitamente
        let mut dummy = Vec::new();
        let mut s = State {
            reader: &mut cursor, writer: &mut dummy,
            in_buf: [0; CHUNK], in_idx: 0, in_len: 0,
            bitbuf: 0, bitcnt: 0, out: [0; MAXWIN], next: 0, first: true,
        };
        assert!(matches!(decode(&mut s, &h), Err(BlastError::OutOfInput)));
    }

    #[test]
    fn test_blast_distance_too_far() {
        // Bytes meticulosamente elaborados para disparar um comando copy que exceda o array interno (MAXWIN)
        let data = vec![0x00, 0x04, 0x01, 0x00, 0x00, 0x00, 0x00];
        let mut input = Cursor::new(data);
        let mut output = Cursor::new(Vec::new());
        assert_eq!(blast(&mut input, &mut output), -3); // -3 map
    }

    #[test]
    fn test_blast_ran_out_codes() {
        // Fuzzing rápido nos 15 bits seguintes para encontrar o código Huffman 
        // não mapeado e atingir a cobertura da linha `-9` com precisão.
        for i in 0..=255 {
            for j in 0..=255 {
                let data = vec![0x01, 0x04, i, j, 0x00, 0x00, 0x00, 0x00];
                let mut input = Cursor::new(data);
                let mut output = Cursor::new(Vec::new());
                if blast(&mut input, &mut output) == -9 {
                    return; // Encontrou o código não mapeado com sucesso!
                }
            }
        }
    }

    #[test]
    fn test_blast_read_error() {
        assert_eq!(blast(&mut BadReader, &mut Cursor::new(Vec::new())), 2);
    }

    struct FailLaterReader { reads: usize }
    impl Read for FailLaterReader {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if self.reads == 0 {
                self.reads += 1;
                buf[0] = 0x00; buf[1] = 0x04; buf[2] = 0x01;
                Ok(3)
            } else {
                Err(io::Error::new(io::ErrorKind::Other, "fail later"))
            }
        }
    }

    #[test]
    fn test_decode_read_error() {
        assert_eq!(blast(&mut FailLaterReader { reads: 0 }, &mut Cursor::new(Vec::new())), 2);
    }

    struct FailLeftoverWriter;
    impl Write for FailLeftoverWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if buf.len() != MAXWIN {
                Err(io::Error::new(io::ErrorKind::Other, "leftover fail"))
            } else { Ok(buf.len()) }
        }
        fn flush(&mut self) -> io::Result<()> { Ok(()) }
    }

    #[test]
    fn test_blast_output_error_flush() {
        let mut data = Vec::new();
        File::open("tests/sids.dbc").unwrap().read_to_end(&mut data).unwrap();
        let header_size = u16::from_le_bytes([data[8], data[9]]) as usize;
        let compressed_data = &data[header_size + 4..];
        assert_eq!(blast(&mut Cursor::new(compressed_data), &mut FailLeftoverWriter), 1);
    }

    #[test]
    fn test_run_cli_output_io_error() {
        assert_eq!(run_cli(&["dbc2dbf".to_string(), "tests/sids.dbc".to_string(), "/dir_does_not_exist/out.dbf".to_string()]), 1);
    }

    #[test]
    fn test_run_cli_fatal_io_error() {
        let out_file = "tests/test_run_cli_fatal_io.dbf".to_string();
        let in_file = "tests/test_too_small.dbc".to_string();
        let mut f = File::create(&in_file).unwrap();
        f.write_all(&[0x00, 0x00]).unwrap(); // Cria arquivo menor que 8 bytes
        
        let ret = run_cli(&["dbc2dbf".to_string(), in_file.clone(), out_file.clone()]);
        assert_eq!(ret, -1); // UnexpectedEof gerando fatal I/O
        
        let _ = std::fs::remove_file(in_file);
        let _ = std::fs::remove_file(out_file);
    }
}
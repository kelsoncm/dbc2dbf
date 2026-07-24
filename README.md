# dbc2dbf 🦀

`dbc2dbf` é um utilitário de linha de comando de alta velocidade escrito em Rust para converter e descompactar arquivos DATASUS DBC no formato de banco de dados DBF.

O formato DBC (utilizado amplamente pelo Ministério da Saúde / DATASUS no Brasil) consiste em arquivos DBF compactados com o algoritmo "implode" da PKWare's Data Compression Library (DCL).

## Uso

```bash
# Compilar via Cargo
cargo build --release

# Executar a conversão de DBC para DBF
./target/release/dbc2dbf input.dbc output.dbf
```

---

## Desenvolvimento Local & Documentação

- **Desenvolvimento Local:** Para orquestração do ambiente de desenvolvimento local, utilize o repositório privado **[workspace](git@github.com:abrasileirado/workspace.git)**.
- **Documentação da Suíte:** Acesse o portal oficial da organização em **[https://abrasileirado.github.io](https://abrasileirado.github.io)**.
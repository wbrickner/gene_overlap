use clap::Clap;
use std::path::PathBuf;

#[derive(Clap, Debug)]
pub struct Arguments {
  #[clap(parse(from_os_str))]
  pub genes: PathBuf,

  #[clap(parse(from_os_str))]
  pub reads: PathBuf,

  #[clap(parse(from_os_str))]
  pub output: PathBuf
}
pub fn parse() -> Arguments { Arguments::parse() }
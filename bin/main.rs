use barbell::search::searcher::Demuxer;

fn main() {
    let mut demuxer = Demuxer::new(0.5);
    let read = b"GGGGGAAATTTGGGCCCCCCCCCCCCCCCCCCCCCC".to_vec();
    let matches = demuxer.demux(&read);
    println!("Matches: {:?}", matches);
}

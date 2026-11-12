//! Read names in the style of Illumina's FASTQ files:
//! `@<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>`.
//!
//! The fragments of a run are laid out like the clusters of one lane of a
//! patterned flow cell as they appear in a demultiplexed FASTQ file: tile
//! numbers increase through the file, within a tile the y coordinate
//! increases in discrete rows, and within a row the x coordinate increases.
//! The name of a fragment depends only on its position in the output, on the
//! total number of fragments and on one random number of the chunk's stream,
//! so it does not depend on the number of threads.

use crate::rng::Rng;

/// Tile numbers of one lane: 2 surfaces × 6 swaths × 78 tiles
/// (1101 … 1178, 1201 … 1678, 2101 … 2678).
pub fn default_tiles() -> Vec<u32> {
    let mut tiles = Vec::with_capacity(2 * 6 * 78);
    for surface in 1..=2u32 {
        for swath in 1..=6u32 {
            for tile in 1..=78u32 {
                tiles.push(surface * 1000 + swath * 100 + tile);
            }
        }
    }
    tiles
}

/// Number of cluster rows per tile.
const ROWS: u64 = 2297;
/// Distance between two rows, in thousandths of a coordinate unit.
const ROW_PITCH_MILLI: u64 = 15_655;
/// First coordinate value in x and y.
const ORIGIN: u64 = 1000;
/// Number of x positions per row.
const X_SPAN: u64 = 31_000;
/// Number of cluster positions per tile, counted row by row.
const POSITIONS: u64 = ROWS * X_SPAN;

/// The largest number of fragments of one run. Up to this number every
/// fragment gets a different name (a tile never holds more fragments than it
/// has positions) and the layout arithmetic cannot overflow.
pub const MAX_FRAGMENTS: u64 = 50_000_000_000;

/// The layout of a whole run.
///
/// Fragment `i` of `total` is on tile `floor(i * tiles / total)`. The
/// positions of a tile (rows of x values, row after row) are divided into as
/// many equal slots as the tile has fragments, and each fragment sits at a
/// random position of its own slot. Coordinates therefore increase through
/// the file, are irregular, and never repeat.
#[derive(Clone, Debug)]
pub struct NameLayout {
    tiles: Vec<u32>,
    total: u64,
}

/// `ceil(a * b / c)`.
#[inline]
fn ceil_mul_div(a: u64, b: u64, c: u64) -> u64 {
    (a * b).div_ceil(c)
}

impl NameLayout {
    pub fn new(total: u64) -> Self {
        Self::with_tiles(total, default_tiles())
    }

    pub fn with_tiles(total: u64, tiles: Vec<u32>) -> Self {
        assert!(!tiles.is_empty());
        assert!((1..=MAX_FRAGMENTS).contains(&total));
        NameLayout { tiles, total }
    }

    /// A cursor positioned on fragment `index` (0-based).
    pub fn cursor(&self, index: u64) -> NameCursor<'_> {
        let mut c = NameCursor {
            layout: self,
            index,
            tile: 0,
            tile_start: 0,
            tile_end: 0,
            slot_start: 0,
        };
        c.enter_tile();
        c
    }
}

/// Walks through the fragments in output order.
#[derive(Clone, Debug)]
pub struct NameCursor<'a> {
    layout: &'a NameLayout,
    index: u64,
    tile: usize,
    /// Fragments `tile_start..tile_end` are on the current tile.
    tile_start: u64,
    tile_end: u64,
    /// First position of the slot of the current fragment.
    slot_start: u64,
}

impl NameCursor<'_> {
    fn enter_tile(&mut self) {
        let n_tiles = self.layout.tiles.len() as u64;
        let total = self.layout.total;
        let k = self.index * n_tiles / total;
        self.tile = k as usize;
        self.tile_start = ceil_mul_div(k, total, n_tiles);
        self.tile_end = ceil_mul_div(k + 1, total, n_tiles);
        let in_tile = self.tile_end - self.tile_start;
        self.slot_start = (self.index - self.tile_start) * POSITIONS / in_tile;
    }

    /// Returns `(tile, x, y)` of the current fragment and moves to the next
    /// one. Uses exactly one random number.
    #[inline]
    pub fn next_coordinates(&mut self, rng: &mut Rng) -> (u32, u32, u32) {
        if self.index >= self.tile_end {
            self.enter_tile();
        }
        let in_tile = self.tile_end - self.tile_start;
        let slot_end = (self.index - self.tile_start + 1) * POSITIONS / in_tile;
        let position = self.slot_start + rng.below((slot_end - self.slot_start).max(1));
        self.slot_start = slot_end;
        self.index += 1;
        let row = position / X_SPAN;
        let x = ORIGIN + position % X_SPAN;
        let y = ORIGIN + (row * ROW_PITCH_MILLI + 500) / 1000;
        (self.layout.tiles[self.tile], x as u32, y as u32)
    }
}

/// Appends the decimal representation of `value` to `out`.
#[inline]
pub fn push_decimal(out: &mut Vec<u8>, mut value: u64) {
    let mut buf = [0u8; 20];
    let mut i = buf.len();
    loop {
        i -= 1;
        buf[i] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }
    out.extend_from_slice(&buf[i..]);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all(total: u64, tiles: Vec<u32>) -> Vec<(u32, u32, u32)> {
        let layout = NameLayout::with_tiles(total, tiles);
        let mut rng = Rng::for_stream(1, 0);
        let mut cursor = layout.cursor(0);
        (0..total)
            .map(|_| cursor.next_coordinates(&mut rng))
            .collect()
    }

    #[test]
    fn tiles_are_as_on_a_two_surface_six_swath_lane() {
        let t = default_tiles();
        assert_eq!(t.len(), 936);
        assert_eq!(t[0], 1101);
        assert_eq!(t[77], 1178);
        assert_eq!(t[78], 1201);
        assert_eq!(*t.last().unwrap(), 2678);
    }

    #[test]
    fn coordinates_increase_and_are_unique() {
        for total in [1u64, 7, 936, 937, 5_000, 300_000] {
            let coords = all(total, default_tiles());
            assert_eq!(coords.len() as u64, total);
            for w in coords.windows(2) {
                let (t0, x0, y0) = w[0];
                let (t1, x1, y1) = w[1];
                assert!(
                    (t0, y0, x0) < (t1, y1, x1),
                    "not increasing: {:?} then {:?} (total {total})",
                    w[0],
                    w[1]
                );
            }
            for &(_, x, y) in &coords {
                assert!((1000..32_000).contains(&x));
                assert!((1000..37_000).contains(&y));
            }
        }
    }

    #[test]
    fn many_fragments_per_row_stay_ordered() {
        // one tile, so that rows hold several fragments
        let coords = all(50_000, vec![1101]);
        let mut sorted = coords.clone();
        sorted.sort_by_key(|&(t, x, y)| (t, y, x));
        sorted.dedup();
        assert_eq!(sorted.len(), coords.len());
        assert_eq!(sorted, coords);
        let max_y = coords.iter().map(|c| c.2).max().unwrap();
        assert!(max_y > 36_000);
    }

    #[test]
    fn all_tiles_are_used_evenly() {
        let coords = all(93_600, default_tiles());
        let mut per_tile = std::collections::BTreeMap::new();
        for (t, _, _) in coords {
            *per_tile.entry(t).or_insert(0u32) += 1;
        }
        assert_eq!(per_tile.len(), 936);
        assert!(per_tile.values().all(|&n| n == 100));
    }

    #[test]
    fn a_cursor_can_start_anywhere() {
        let total = 20_000u64;
        let layout = NameLayout::new(total);
        // With the same random numbers, a cursor started in the middle
        // gives the same coordinates as one that walked there.
        let draws: Vec<Rng> = (0..total).map(|i| Rng::for_stream(99, i)).collect();
        let mut from_start = layout.cursor(0);
        let reference: Vec<(u32, u32, u32)> = (0..total)
            .map(|i| from_start.next_coordinates(&mut draws[i as usize].clone()))
            .collect();
        for start in [0u64, 1, 21, 22, 4096, 8191, 19_999] {
            let mut cursor = layout.cursor(start);
            for i in start..(start + 50).min(total) {
                let c = cursor.next_coordinates(&mut draws[i as usize].clone());
                assert_eq!(c, reference[i as usize], "fragment {i} from {start}");
            }
        }
    }

    #[test]
    fn decimal_formatting() {
        let mut v = Vec::new();
        push_decimal(&mut v, 0);
        v.push(b' ');
        push_decimal(&mut v, 1101);
        v.push(b' ');
        push_decimal(&mut v, u64::MAX);
        assert_eq!(v, b"0 1101 18446744073709551615");
    }
}

//! Simple in-memory key-value cache.
//!
//! Wraps a [`std::collections::HashMap`] with a stable API so the backing
//! store can be swapped for `redb` or `sled` without touching call sites.

use std::collections::HashMap;

/// In-memory key-value store that maps `String` keys to raw byte values.
///
/// Values are stored as `Vec<u8>` so that any serialised type can be cached
/// without coupling to a specific codec.
pub struct KvStore {
    inner: HashMap<String, Vec<u8>>,
}

impl KvStore {
    /// Creates an empty [`KvStore`].
    pub fn new() -> Self {
        Self {
            inner: HashMap::new(),
        }
    }

    /// Returns a reference to the value stored under `key`, or `None`.
    pub fn get(&self, key: &str) -> Option<&[u8]> {
        self.inner.get(key).map(Vec::as_slice)
    }

    /// Inserts or replaces the value for `key`.
    pub fn set(&mut self, key: String, value: Vec<u8>) {
        self.inner.insert(key, value);
    }

    /// Removes and returns the value for `key`, or `None` if absent.
    pub fn remove(&mut self, key: &str) -> Option<Vec<u8>> {
        self.inner.remove(key)
    }
}

impl Default for KvStore {
    fn default() -> Self {
        Self::new()
    }
}

/**
 * ProteinScope Mobile Service Worker
 * Provides 100% offline caching for all mobile assets, local 3Dmol library, and 9 cached proteins.
 */

const CACHE_NAME = 'proteinscope-mobile-v2';
const ASSETS_TO_CACHE = [
    './',
    './index.html',
    './manifest.json',
    './icons/icon.svg',
    './js/3Dmol-min.js',
    './js/cached_data.js',
    './js/biophysics_engine.js',
    './js/app.js'
];

self.addEventListener('install', (e) => {
    e.waitUntil(
        caches.open(CACHE_NAME).then((cache) => {
            console.log('[ServiceWorker] Pre-caching offline workstation assets');
            return cache.addAll(ASSETS_TO_CACHE).catch(err => {
                console.warn('[ServiceWorker] Some assets could not be cached immediately:', err);
            });
        })
    );
    self.skipWaiting();
});

self.addEventListener('activate', (e) => {
    e.waitUntil(
        caches.keys().then((keys) => {
            return Promise.all(
                keys.map((key) => {
                    if (key !== CACHE_NAME) {
                        console.log('[ServiceWorker] Purging old cache:', key);
                        return caches.delete(key);
                    }
                })
            );
        })
    );
    self.clients.claim();
});

self.addEventListener('fetch', (e) => {
    // For online API requests (PDB, UniProt, NCBI, AlphaFold), try network first
    if (
        e.request.url.includes('ncbi.nlm.nih.gov') || 
        e.request.url.includes('alphafold.ebi.ac.uk') ||
        e.request.url.includes('rcsb.org') ||
        e.request.url.includes('uniprot.org')
    ) {
        e.respondWith(
            fetch(e.request).catch(() => {
                return new Response(JSON.stringify({ error: 'Device is currently offline. Online retrieval requires an internet connection.' }), {
                    headers: { 'Content-Type': 'application/json' }
                });
            })
        );
        return;
    }

    // For static app assets, cache first, fallback to network
    e.respondWith(
        caches.match(e.request).then((cachedResponse) => {
            if (cachedResponse) {
                return cachedResponse;
            }
            return fetch(e.request).then((networkResponse) => {
                if (networkResponse && networkResponse.status === 200 && networkResponse.type === 'basic') {
                    const responseClone = networkResponse.clone();
                    caches.open(CACHE_NAME).then((cache) => {
                        cache.put(e.request, responseClone);
                    });
                }
                return networkResponse;
            });
        })
    );
});

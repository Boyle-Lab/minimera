;;; Copyright 2025 Steve Losh and contributors
;;; SPDX-License-Identifier: GPL-3.0-or-later

(in-package :minimera)

;;;; Hashing ------------------------------------------------------------------

;;; https://naml.us/post/inverse-of-a-hash-function/

(deftype u64 ()
  `(integer 0 ,(1- (expt 2 64))))

(defun-inline wrap64 (n)
  (mod n (expt 2 64)))

(defun-inline shl (x bits)
  (logand (ash x bits)
          (1- (ash 1 64))))

(defun-inline shr (x bits)
  (logand (ash x (- bits))
          (1- (ash 1 64))))

(defun-inline hash64 (n)
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (etypecase n
    (fixnum
      (_ n
        (wrap64 (+ (lognot _) (shl _ 21)))
        (wrap64 (logxor _ (shr _ 24)))
        (wrap64 (+ (+ _ (shl _ 3))
                   (shl _ 8)))
        (wrap64 (logxor _ (shr _ 14)))
        (wrap64 (+ (+ _ (shl _ 2))
                   (shl _ 4)))
        (wrap64 (logxor _ (shr _ 28)))
        (wrap64 (+ _ (shl _ 31)))))))

(defun ihash64 (n)
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  ;;; https://naml.us/post/inverse-of-a-hash-function/
  (etypecase n
    (fixnum
      (let* ((key n)

             ;;   // Invert key = key + (key << 31)
             ;;   tmp = key - (key <<31);
             ;;   key = key - (tmp <<31);
             (tmp (wrap64 (- key (shl key 31))))
             (key (wrap64 (- key (shl tmp 31))))

             ;;   // Invert key = key ^ (key >> 28)
             ;;   tmp = key ^ key >>28;
             ;;   key = key ^ tmp >>28;
             (tmp (wrap64 (logxor key (shr key 28))))
             (key (wrap64 (logxor key (shr tmp 28))))

             ;;   // Invert key *= 21
             ;;   key *= 14933078535860113213u;
             (key (wrap64 (* key 14933078535860113213)))

             ;;   // Invert key = key ^ (key >> 14)
             ;;   tmp = key ^ key >> 14;
             ;;   tmp = key ^ tmp >> 14;
             ;;   tmp = key ^ tmp >> 14;
             ;;   key = key ^ tmp >> 14;
             (tmp (wrap64 (logxor key (shr key 14))))
             (tmp (wrap64 (logxor key (shr tmp 14))))
             (tmp (wrap64 (logxor key (shr tmp 14))))
             (key (wrap64 (logxor key (shr tmp 14))))

             ;;   // Invert key *= 265
             ;;   key *= 15244667743933553977u;
             (key (wrap64 (* key 15244667743933553977)))

             ;;   // Invert key = key ^ (key >> 24)
             ;;   tmp = key ^ key >> 24;
             ;;   key = key ^ tmp >> 24;
             (tmp (wrap64 (logxor key (shr key 24))))
             (key (wrap64 (logxor key (shr tmp 24))))

             ;;   // Invert key = (~key) + (key << 21)
             ;;   tmp = ~key;
             ;;   tmp = ~( key -( tmp <<21));
             ;;   tmp = ~( key -( tmp <<21));
             ;;   key = ~( key -( tmp <<21));
             (tmp (wrap64 (lognot key)))
             (tmp (wrap64 (lognot (- key (shl tmp 21)))))
             (tmp (wrap64 (lognot (- key (shl tmp 21)))))
             (key (wrap64 (lognot (- key (shl tmp 21))))))
        key))))

(defun phihash-string (string &key (start 0) (end (length string)))
  "Return the φ-hash of the subsequence of `string` bounded by `start` and `end`."
  (iterate
    (with hash = 0)
    (for ch :in-string string :from start :below end)
    (zapf hash (_ %
                 (ash _ 2)
                 (logior _ (ecase ch
                             (#\A #b00)
                             (#\C #b01)
                             (#\G #b10)
                             (#\T #b11)))))
    (returning (hash64 hash))))

(defun-inline phihash-chunk (chunk)
  (hash64 chunk))


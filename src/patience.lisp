(in-package :minimera)

;;;; Longest Increasing Subsequence -------------------------------------------

;;; Based on Rosetta code, cleaned up and extended to allow :key function:
;;;
;;; https://rosettacode.org/wiki/Longest_increasing_subsequence#Common_Lisp:_Using_the_Patience_Sort_approach

(defun patience-insert (item piles &key key)
  ;; piles is an array that will be mutated.  each pile is ((item . backpointer)
  ;; …) where backpointer points at the previous element (and this conveniently
  ;; happens to resemble vanilla list structure).
  (loop :for pile :across piles
        :and prev = nil :then pile
        :for i from 0
        :for top = (car pile)
        :do (when (<= (funcall key item)
                      (funcall key (car top)))
              (push (cons item (car prev)) (aref piles i))
              (return))
        :finally (vector-push-extend (list (cons item top)) piles)))

(defun longest-increasing-subsequence (sequence &key (key #'identity) (result-type 'list))
  "Return the longest increasing subsequence of elements from `sequence`.

  Elements will be compared with `<`.  `key` will be called on each first.

  The result will be a sequence of type `result-type`.

  "
  (let ((piles (make-array 256 :adjustable t :fill-pointer 0)))
    (map nil (lambda (item) (patience-insert item piles :key key)) sequence)
    (nreverse (coerce (car (aref piles (1- (length piles)))) result-type))))

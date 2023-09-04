#define       N              %r15
#define       N2             %r14
#define       N1             %r13
#define       X              %r12
#define       W              %r11
#define       k2             %r10
#define       I              %r9
#define       NLO            %r8
#define       w0             %rbx
#define       level          -8(%rbp)
#define       NTW            -16(%rbp)
#define       NLO2           -24(%rbp)
#define       NLO4           -32(%rbp)
#define       sptr           -40(%rbp)
#define       STACK_SIZE     $64

              .global        fft_1d_batch_dif

              .text

fft_1d_batch_dif:
              push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %rbx
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              mov            %rdi, X
			  mov            %rsi, W
			  mov            %rdx, N1
              mov            %rcx, N
              mov            %r8, NLO
              xor            %rax, %rax
              mov            %rax, level
              mov            $1, %rax
              mov            %rax, NTW
              mov            NLO, %rax
              mov            NLO, %rdx
              and            $0xfffffffffffffffc, %rax
              and            $0xfffffffffffffffe, %rdx
              mov            %rax, NLO4
              mov            %rdx, NLO2
              imul           $64, N1, %rax
              mov            %rsp, sptr
              sub            %rax, %rsp
              mov            %rsp, w0
              call           next_level
              mov            sptr, %rsp
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %rbx
              mov            %rbp, %rsp
              pop            %rbp
              ret
next_level:   xor            %rdx, %rdx
              mov            N, %rax
              divq           N1
              mov            %rax, N2
              xor            k2, k2
k2loop:       mov            $1, I
w4:           mov            NTW, %rax
              imul           k2, %rax
              imul           I, %rax
              shl            %rax
              vbroadcastsd   (W, %rax, 8), %ymm0
              vbroadcastsd   8(W, %rax, 8), %ymm1
              mov            I, %rax
              shl            $3, %rax
              vmovupd        %ymm0, (w0, %rax, 8)
              vmovupd        %ymm1, 32(w0, %rax, 8)
              inc            I
              cmp            I, N1
              jne            w4
              xor            I, I
iloop_begin:  cmp            I, NLO4
              jle            iloop_end
              mov            k2, %rax
              imul           NLO, %rax
              add            I, %rax
              shl            %rax
              cmp            $0, level
              jne            skip_shuf4
              xor            %rdx, %rdx
shuf4:        mov            %rdx, %rcx
              imul           N2, %rcx
              imul           NLO, %rcx
              shl            %rcx
              add            %rax, %rcx
              vpermpd        $216, (X, %rcx, 8), %ymm0
              vpermpd        $216, 32(X, %rcx, 8), %ymm1
              vperm2f128     $0x20, %ymm1, %ymm0, %ymm2
              vperm2f128     $0x31, %ymm1, %ymm0, %ymm3
              vmovupd        %ymm2, (X, %rcx, 8)
              vmovupd        %ymm3, 32(X, %rcx, 8)
              inc            %rdx
              cmp            %rdx, N1
              jne            shuf4
skip_shuf4:   lea            (X, %rax, 8), %rdi
              mov            %rdi, %rsi
              add            $32, %rsi
              mov            NLO, %rdx
              imul           N2, %rdx
              shl            %rdx
              mov            N1, %rcx
              push           %r8
              push           %r9
              push           %r10
              push           %r11
              mov            w0, %r8
              call           sfft_complex_dif_w4
              pop            %r11
              pop            %r10
              pop            %r9
              pop            %r8
	          add            $4, I
              jmp            iloop_begin
iloop_end:    mov            NLO4, %rax
              cmp            %rax, NLO2
              je             skip2
              mov            $1, I
w2:           mov            NTW, %rax
              imul           k2, %rax
              imul           I, %rax
              shl            %rax
              vpbroadcastq   (W, %rax, 8), %xmm0
              vpbroadcastq   8(W, %rax, 8), %xmm1
              mov            I, %rax
              shl            $2, %rax
              vmovupd        %xmm0, (w0, %rax, 8)
              vmovupd        %xmm1, 16(w0, %rax, 8)
              inc            I
              cmp            I, N1
              jne            w2
              mov            k2, %rax
              imul           NLO, %rax
              add            NLO4, %rax
              shl            %rax
              cmp            $0, level
              jne            skip_shuf2
              xor            %rdx, %rdx
shuf2:        mov            %rdx, %rcx
              imul           N2, %rcx
              imul           NLO, %rcx
              shl            %rcx
              add            %rax, %rcx
              vmovq          8(X, %rcx, 8), %xmm0
              vmovq          16(X, %rcx, 8), %xmm1
              vmovq          %xmm1, 8(X, %rcx, 8)
              vmovq          %xmm0, 16(X, %rcx, 8)
              inc            %rdx
              cmp            %rdx, N1
              jne            shuf2
skip_shuf2:   lea            (X, %rax, 8), %rdi
              mov            %rdi, %rsi
              add            $16, %rsi
              mov            NLO, %rdx
              imul           N2, %rdx
              shl            %rdx
              mov            N1, %rcx
              push           %r8
              push           %r9
              push           %r10
              push           %r11
              mov            w0, %r8
              call           sfft_complex_dif_w2
              pop            %r11
              pop            %r10
              pop            %r9
              pop            %r8
skip2:        mov            NLO2, %rax
              cmp            %rax, NLO
              je             skip1
              mov            $1, I
w1:           mov            NTW, %rax
              imul           I, %rax
              imul           k2, %rax
              shl            %rax
              mov            (W, %rax, 8), %rcx
              mov            8(W, %rax, 8), %rdx
              mov            I, %rax
              shl            %rax
              mov            %rcx, (w0, %rax, 8)
              mov            %rdx, 8(w0, %rax, 8)
              inc            I
              cmp            I, N1
              jne            w1
              mov            k2, %rax
              imul           NLO, %rax
              add            NLO2, %rax
              shl            %rax
              lea            (X, %rax, 8), %rdi
              mov            %rdi, %rsi
              add            $8, %rsi
              mov            NLO, %rdx
              imul           N2, %rdx
              shl            %rdx
              mov            N1, %rcx
              push           %r8
              push           %r9
              push           %r10
              push           %r11
              mov            w0, %r8
              call           sfft_complex_dif_w1
              pop            %r11
              pop            %r10
              pop            %r9
              pop            %r8
skip1:        inc            k2
              cmp            k2, N2
              jne            k2loop
              mov            N1, %rax
              cmp            %rax, N
              je             done
              xor            I, I
next_loop:    push           X
              push           NTW
              push           N
              push           N2
              push           I
              push           level
              mov            N2, %rax
              imul           NLO, %rax
              imul           I, %rax
              shl            %rax
              lea            (X, %rax, 8), X
              mov            NTW, %rax
              imul           N1, %rax
              mov            %rax, NTW
              mov            N2, %rdx
              mov            %rdx, N
              mov            level, %rax
              inc            %rax
              mov            %rax, level
              call           next_level
              pop            level
              pop            I
              pop            N2
              pop            N
              pop            NTW
              pop            X
              inc            I
              cmp            I, N1
              jne            next_loop
done:         ret

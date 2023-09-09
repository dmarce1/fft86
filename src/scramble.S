#define       X              %r8
#define       R              %r9
#define       N              %r10
#define       M              %r11
#define       Rm1            %r12
#define       Nm1            %r13
#define       NoR            %r14
#define       I              %r15
#define       J              %rdi
#define       K              %rsi
#define       m              %rcx
#define       M2             -8(%rbp)
#define       M4             -16(%rbp)
#define       M32            -24(%rbp)
#define       STACK_SIZE     $32

              .global        scramble_hi

              .text


scramble_hi:  push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              mov            %rdi, X
              mov            %rsi, R
              mov            %rdx, N
              mov            %rcx, M
              mov            M, %rax
              mov            M, %rdx
              mov            M, %rcx
              and            $0xfffffffffffffffc, %rax
              and            $0xfffffffffffffffe, %rdx
              and            $0xffffffffffffffc0, %rcx
              mov            %rax, M4
              mov            %rdx, M2
              mov            %rcx, M32
              xor            %edx, %edx
              mov            N, %rax
              div            %r9d
              mov            %rax, NoR
              mov            N, Nm1
              dec            Nm1
              mov            R, Rm1
              dec            Rm1
              xor            I, I
              xor            J, J
mainloop:     cmp            I, J
              jge            skip
              xor            m, m
              mov            I, %rax
              mov            J, %rdx
              imul           M, %rax
              imul           M, %rdx
              cmp            $0, M32
              je             small_begin
big_mloop:    vmovupd        (X, %rax, 8), %ymm0
              vmovupd        32(X, %rax, 8), %ymm1
              vmovupd        64(X, %rax, 8), %ymm2
              vmovupd        96(X, %rax, 8), %ymm3
              vmovupd        128(X, %rax, 8), %ymm4
              vmovupd        160(X, %rax, 8), %ymm5
              vmovupd        192(X, %rax, 8), %ymm6
              vmovupd        224(X, %rax, 8), %ymm7
              vmovupd        (X, %rdx, 8), %ymm8
              vmovupd        32(X, %rdx, 8), %ymm9
              vmovupd        64(X, %rdx, 8), %ymm10
              vmovupd        96(X, %rdx, 8), %ymm11
              vmovupd        128(X, %rdx, 8), %ymm12
              vmovupd        160(X, %rdx, 8), %ymm13
              vmovupd        192(X, %rdx, 8), %ymm14
              vmovupd        224(X, %rdx, 8), %ymm15
              vmovupd        %ymm0, (X, %rdx, 8)
              vmovupd        %ymm1, 32(X, %rdx, 8)
              vmovupd        %ymm2, 64(X, %rdx, 8)
              vmovupd        %ymm3, 96(X, %rdx, 8)
              vmovupd        %ymm4, 128(X, %rdx, 8)
              vmovupd        %ymm5, 160(X, %rdx, 8)
              vmovupd        %ymm6, 192(X, %rdx, 8)
              vmovupd        %ymm7, 224(X, %rdx, 8)
              vmovupd        %ymm8, (X, %rax, 8)
              vmovupd        %ymm9, 32(X, %rax, 8)
              vmovupd        %ymm10, 64(X, %rax, 8)
              vmovupd        %ymm11, 96(X, %rax, 8)
              vmovupd        %ymm12, 128(X, %rax, 8)
              vmovupd        %ymm13, 160(X, %rax, 8)
              vmovupd        %ymm14, 192(X, %rax, 8)
              vmovupd        %ymm15, 224(X, %rax, 8)
              add            $32, %rax
              add            $32, %rdx
              add            $32, m
              cmp            m, M32
              jne            big_mloop
small_begin:  cmp            m, M4
              jle            test2
              vmovupd        (X, %rax, 8), %ymm0
              vmovupd        (X, %rdx, 8), %ymm1
              vmovupd        %ymm1, (X, %rax, 8)
              vmovupd        %ymm0, (X, %rdx, 8)
              add            $4, %rax
              add            $4, %rdx
              add            $4, m
              jmp            small_begin
test2:        mov            M2, %rax
              cmp            %rax, M4
              je             no2
              mov            I, %rax
              mov            J, %rdx
              imul           M, %rax
              imul           M, %rdx
              add            M4, %rax
              add            M4, %rdx
              vmovupd        (X, %rax, 8), %xmm0
              vmovupd        (X, %rdx, 8), %xmm1
              vmovupd        %xmm1, (X, %rax, 8)
              vmovupd        %xmm0, (X, %rdx, 8)
no2:          cmp            M, M2
              je             skip
              mov            I, %rax
              mov            J, %rdx
              imul           M, %rax
              imul           M, %rdx
              add            M2, %rax
              add            M2, %rdx
              vmovq          (X, %rax, 8), %xmm0
              vmovq          (X, %rdx, 8), %xmm1
              vmovq          %xmm1, (X, %rax, 8)
              vmovq          %xmm0, (X, %rdx, 8)
skip:         cmp            I, Nm1
              je             done
              mov            NoR, K
              inc            J
rev_loop:     mov            K, %rax
              imul           Rm1, %rax
              cmp            %rax, J
              jle            rev_done
              sub            %rax, J
              xor            %edx, %edx
              mov            K, %rax
              div            %r9d
              mov            %rax, K
              jmp            rev_loop
rev_done:     add            K, J
              dec            J
              inc            I
              jmp            mainloop
done:         pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              mov            %rbp, %rsp
              pop            %rbp
              ret

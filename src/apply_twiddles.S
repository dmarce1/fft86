#define       X              %r15
#define       Wr             %r14
#define       Wi             %r13
#define       R              %r12
#define       N1             %r11
#define       N2             %r10
#define       k2             %r9
#define       n1             %r8
#define       k2rev          %rdi
#define       I              %rsi
#define       N2oR           -8(%rbp)
#define       Rm1            -16(%rbp)
#define       Nm1            -24(%rbp)
#define       N12            -32(%rbp)
#define       N14            -40(%rbp)
#define       STACK_SIZE     $128



              .global        apply_twiddles

              .text
apply_twiddles:
              push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %r12
              push           %r13
              push           %r14
              push           %r15
              mov            %rdi, X
              mov            %rsi, Wr
              mov            %rdx, Wi
              mov            %rcx, R
              mov            %r8, N1
              mov            %r9, N2
              xor            %rdx, %rdx
              mov            N2, %rax
              div            R
              mov            %rax, N2oR
              mov            R, %rdx
              dec            %rdx
              mov            %rdx, Rm1
              mov            N2, %rax
              dec            %rax
              mov            %rax, Nm1
              mov            N1, %rax
              mov            N1, %rdx
              and            $0xfffffffffffffffc, %rax
              and            $0xfffffffffffffffe, %rdx
              mov            %rax, N14
              mov            %rdx, N12
              xor            k2, k2
              xor            k2rev, k2rev
k2loop:       xor            n1, n1
n1loop_begin: cmp            n1, N14
              jle            n1loop_end
              mov            k2, %rdx
              imul           N1, %rdx
              add            n1, %rdx
              mov            k2rev, %rax
              imul           N1, %rax
              add            n1, %rax
              shl            %rax
              vmovupd        (Wr, %rdx, 8), %ymm10
              vmovupd        (Wi, %rdx, 8), %ymm11
              vmovupd        (X, %rax, 8), %ymm0
              vmovupd        32(X, %rax, 8), %ymm1
              vmulpd         %ymm11, %ymm1, %ymm2
              vmulpd         %ymm10, %ymm1, %ymm3
              vfmsub231pd    %ymm10, %ymm0, %ymm2
              vfmadd231pd    %ymm11, %ymm0, %ymm3
              vperm2f128     $0x20, %ymm3, %ymm2, %ymm0
              vperm2f128     $0x31, %ymm3, %ymm2, %ymm1
              vpermpd        $216, %ymm0, %ymm0
              vpermpd        $216, %ymm1, %ymm1
              vmovupd        %ymm0, (X, %rax, 8)
              vmovupd        %ymm1, 32(X, %rax, 8)
              add            $4, n1
              jmp            n1loop_begin
n1loop_end:   mov            N12, %rax
              cmp            %rax, N14
              je             skip2
              mov            k2, %rdx
              imul           N1, %rdx
              add            N14, %rdx
              mov            k2rev, %rax
              imul           N1, %rax
              add            N14, %rax
              shl            %rax
              vmovupd        (Wr, %rdx, 8), %xmm10
              vmovupd        (Wi, %rdx, 8), %xmm11
              vmovupd        (X, %rax, 8), %xmm0
              vmovupd        16(X, %rax, 8), %xmm1
              vmulpd         %xmm11, %xmm1, %xmm2
              vmulpd         %xmm10, %xmm1, %xmm3
              vfmsub231pd    %xmm10, %xmm0, %xmm2
              vfmadd231pd    %xmm11, %xmm0, %xmm3
              vmovupd        %xmm2, (X, %rax, 8)
              vmovupd        %xmm3, 16(X, %rax, 8)
              vmovq          8(X, %rax, 8), %xmm0
              vmovq          16(X, %rax, 8), %xmm1
              vmovq          %xmm1, 8(X, %rax, 8)
              vmovq          %xmm0, 16(X, %rax, 8)
skip2:        cmp            N1, N12
              je             skip1
              mov            k2, %rdx
              imul           N1, %rdx
              add            N12, %rdx
              mov            k2rev, %rax
              imul           N1, %rax
              add            N12, %rax
              shl            %rax
              vmovq          (Wr, %rdx, 8), %xmm10
              vmovq          (Wi, %rdx, 8), %xmm11
              vmovq          (X, %rax, 8), %xmm0
              vmovq          8(X, %rax, 8), %xmm1
              vmulsd         %xmm11, %xmm1, %xmm2
              vmulsd         %xmm10, %xmm1, %xmm3
              vfmsub231sd    %xmm10, %xmm0, %xmm2
              vfmadd231sd    %xmm11, %xmm0, %xmm3
              vmovq          %xmm2, (X, %rax, 8)
              vmovq          %xmm3, 8(X, %rax, 8)
skip1:        cmp            k2, Nm1
              je             skip_revcntr
              mov            N2oR, I
              inc            k2rev
revcntr_loop: mov            I, %rax
              imul           Rm1, %rax
              cmp            %rax, k2rev
              jle            revcntr_done
              sub            %rax, k2rev
              xor            %rdx, %rdx
              mov            I, %rax
              divq           R
              mov            %rax, I
              jmp            revcntr_loop
revcntr_done: add            I, k2rev
              dec            k2rev
skip_revcntr: inc            k2
              cmp            N2, k2
              jne            k2loop
              pop            %r15
              pop            %r14
              pop            %r13
              pop            %r12
              leave
              ret
              .align         32
IOTA:         .quad          0
              .quad          2
              .quad          4
              .quad          6
ALL_TRUE:     .double        -1
              .double        -1
              .double        -1
              .double        -1
ONE:          .double        1.0
ZERO:         .double        0.0

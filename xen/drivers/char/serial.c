/******************************************************************************
 * serial.c
 * 
 * Framework for serial device drivers.
 * 
 * Copyright (c) 2003-2005, K A Fraser
 */

#include <xen/config.h>
#include <xen/init.h>
#include <xen/irq.h>
#include <xen/keyhandler.h> 
#include <xen/reboot.h>
#include <xen/sched.h>
#include <xen/serial.h>

static struct serial_port com[2] = {
    { .lock = SPIN_LOCK_UNLOCKED }, 
    { .lock = SPIN_LOCK_UNLOCKED }
};

void serial_rx_interrupt(struct serial_port *port, struct cpu_user_regs *regs)
{
    char c;
    serial_rx_fn fn;
    unsigned long flags;

    BUG_ON(!port->driver);
    BUG_ON(!port->driver->getc);

    for ( ; ; )
    {
        spin_lock_irqsave(&port->lock, flags);

        if ( !port->driver->getc(port, &c) )
            break;

        fn = NULL;
        if ( port->rx != NULL )
            fn = port->rx;
        else if ( (c & 0x80) && (port->rx_hi != NULL) )
            fn = port->rx_hi;
        else if ( !(c & 0x80) && (port->rx_lo != NULL) )
            fn = port->rx_lo;
        else if ( (port->rxbufp - port->rxbufc) != RXBUFSZ )
            port->rxbuf[MASK_RXBUF_IDX(port->rxbufp++)] = c;            

        spin_unlock_irqrestore(&port->lock, flags);

        if ( fn != NULL )
            (*fn)(c & 0x7f, regs);

        cpu_relax();
    }

    spin_unlock_irqrestore(&port->lock, flags);
}

void serial_putc(int handle, char c)
{
    struct serial_port *port = &com[handle & SERHND_IDX];
    unsigned long flags;

    if ( (handle == -1) || !port->driver || !port->driver->putc )
        return;

    spin_lock_irqsave(&port->lock, flags);

    if ( (c == '\n') && (handle & SERHND_COOKED) )
        port->driver->putc(port, '\r');

    if ( handle & SERHND_HI )
        c |= 0x80;
    else if ( handle & SERHND_LO )
        c &= 0x7f;

    port->driver->putc(port, c);

    spin_unlock_irqrestore(&port->lock, flags);
}

void serial_puts(int handle, const char *s)
{
    while ( *s != '\0' )
        serial_putc(handle, *s++);
}

char serial_getc(int handle)
{
    struct serial_port *port = &com[handle & SERHND_IDX];
    char c;
    unsigned long flags;

    if ( (handle == -1) || !port->driver || !port->driver->getc )
        return '\0';

    do {        
        for ( ; ; )
        {
            spin_lock_irqsave(&port->lock, flags);
            
            if ( port->rxbufp != port->rxbufc )
            {
                c = port->rxbuf[MASK_RXBUF_IDX(port->rxbufc++)];
                break;
            }
            
            if ( port->driver->getc(port, &c) )
                break;

            spin_unlock_irqrestore(&port->lock, flags);

            cpu_relax();
        }
    } while ( ((handle & SERHND_LO) &&  (c & 0x80)) ||
              ((handle & SERHND_HI) && !(c & 0x80)) );
    
    return c & 0x7f;
}

int serial_parse_handle(char *conf)
{
    int handle;

    /* Silently fail if user has explicitly requested no serial I/O. */
    if ( strcmp(conf, "none") == 0 )
        return -1;

    if ( strncmp(conf, "com", 3) != 0 )
        goto fail;

    switch ( conf[3] )
    {
    case '1':
        handle = 0;
        break;
    case '2':
        handle = 1;
        break;
    default:
        goto fail;
    }

    if ( conf[4] == 'H' )
        handle |= SERHND_HI;
    else if ( conf[4] == 'L' )
        handle |= SERHND_LO;

    handle |= SERHND_COOKED;

    return handle;

 fail:
    printk("ERROR: bad serial-interface specification '%s'\n", conf);
    return -1;
}

void serial_set_rx_handler(int handle, serial_rx_fn fn)
{
    struct serial_port *port = &com[handle & SERHND_IDX];
    unsigned long flags;

    if ( handle == -1 )
        return;

    spin_lock_irqsave(&port->lock, flags);

    if ( port->rx != NULL )
        goto fail;

    if ( handle & SERHND_LO )
    {
        if ( port->rx_lo != NULL )
            goto fail;
        port->rx_lo = fn;        
    }
    else if ( handle & SERHND_HI )
    {
        if ( port->rx_hi != NULL )
            goto fail;
        port->rx_hi = fn;
    }
    else
    {
        if ( (port->rx_hi != NULL) || (port->rx_lo != NULL) )
            goto fail;
        port->rx = fn;
    }

    spin_unlock_irqrestore(&port->lock, flags);
    return;

 fail:
    spin_unlock_irqrestore(&port->lock, flags);
    printk("ERROR: Conflicting receive handlers for COM%d\n", 
           handle & SERHND_IDX);
}

void serial_force_unlock(int handle)
{
    struct serial_port *port = &com[handle & SERHND_IDX];
    if ( handle != -1 )
        port->lock = SPIN_LOCK_UNLOCKED;
}

void serial_init_preirq(void)
{
    int i;
    for ( i = 0; i < ARRAY_SIZE(com); i++ )
        if ( com[i].driver && com[i].driver->init_preirq )
            com[i].driver->init_preirq(&com[i]);
}

void serial_init_postirq(void)
{
    int i;
    for ( i = 0; i < ARRAY_SIZE(com); i++ )
        if ( com[i].driver && com[i].driver->init_postirq )
            com[i].driver->init_postirq(&com[i]);
}

void serial_endboot(void)
{
    int i;
    for ( i = 0; i < ARRAY_SIZE(com); i++ )
        if ( com[i].driver && com[i].driver->endboot )
            com[i].driver->endboot(&com[i]);
}

void serial_register_uart(int idx, struct uart_driver *driver, void *uart)
{
    com[idx].driver = driver;
    com[idx].uart   = uart;
}

/*
 * Local variables:
 * mode: C
 * c-set-style: "BSD"
 * c-basic-offset: 4
 * tab-width: 4
 * indent-tabs-mode: nil
 * End:
 */

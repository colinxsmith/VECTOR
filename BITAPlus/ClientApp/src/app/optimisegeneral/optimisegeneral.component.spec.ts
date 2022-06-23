import { ComponentFixture, TestBed } from '@angular/core/testing';

import { OptimisegeneralComponent } from './optimisegeneral.component';

describe('OptimisegeneralComponent', () => {
  let component: OptimisegeneralComponent;
  let fixture: ComponentFixture<OptimisegeneralComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ OptimisegeneralComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(OptimisegeneralComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});

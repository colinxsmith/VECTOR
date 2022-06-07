import { ComponentFixture, TestBed } from '@angular/core/testing';

import { OptimiseComponent } from './optimise.component';

describe('OptimiseComponent', () => {
  let component: OptimiseComponent;
  let fixture: ComponentFixture<OptimiseComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ OptimiseComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(OptimiseComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
